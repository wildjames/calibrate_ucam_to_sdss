import argparse
import copy
import json
from os import mkdir
from os.path import isdir, join
from pprint import pprint
from time import sleep

import hipercam as hcam
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import coordinates as coord
from astropy import time
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from ruamel import yaml

parser = argparse.ArgumentParser("YAML input method.")
parser.add_argument(
    "input_file",
    help="The input YAML file to be computed"
)
args = parser.parse_args()
yaml_fname = args.input_file

with open(yaml_fname) as yaml_file:
    input_dict = yaml.safe_load(yaml_file)

fnames =                    input_dict['fnames']
oname =                     input_dict['oname']
lc_dir =                    input_dict['lc_dir']
target_coords =             input_dict['target_coords']
T0 =                        input_dict['T0']
period =                    input_dict['period']
lower_phase, upper_phase =  input_dict['phase limits']
k_ext =                     input_dict['k_ext']
obsname =                   input_dict['obsname']
bands =                     input_dict['bands']
comparison_aperture =       input_dict['comparison_aperture']
comparison_mags_tablename = input_dict['comparison_mags_tablename']
variables_fname =           input_dict['previous values']


def tcorrect(tseries, star, observatory, type='B'):
    """
    Correct for light travel time.

    Arguments:
    ----------
    tseries: hipercam.hlog.Tseries
        Time series object

    star: astropy.coordinate.SkyCoord
        Location of star on Sky

    observatory: string
        Observatory name. See coord.EarthLocation.get_site_names() for list. If not in the list, assumed
        to be "lat, lon", comma separated.

    type: string (default=B)
        Heliocentric (H) or Barcentric (B)

    Returns
    -------
    tseries_corr : hipercam.hlog.Tseries
        Time series object with corrected time axis
    """
    ts = copy.deepcopy(tseries)
    try:
        location = coord.EarthLocation.of_site(observatory)
    except:
        lat, lon = observatory.split(',')
        print("  Attempting to get the earth location from latitude and longitude")
        location = coord.EarthLocation.from_geodetic(lat=lat, lon=lon)

    times = time.Time(tseries.t, format='mjd', scale='utc', location=location)

    if type == 'B':
        corr = times.light_travel_time(star, 'barycentric', )
        corr = times.tdb + corr
    else:
        corr = times.light_travel_time(star, 'heliocentric')
        corr = times.utc + corr
    ts.t = corr.mjd
    return ts

def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10.0**(-mag/2.5)
    flux*= alpha

    return flux

##### INFORMATION GATHERING #####

# Get back previously calculated values
with open(variables_fname, 'r') as f:
    variables = json.load(f)
# Colour terms. a_u is u-g, g and r are g-r.
a_u = variables['au']
a_g = variables['ag']
a_r = variables['ar']

# Zero points in each CCD for UCAM on the NTT.
u_zp = variables['u_zp']
g_zp = variables['g_zp']
r_zp = variables['r_zp']


# Gather the aperture magnitudes. Despite the name, this DOES contain the target aperture, too!
comparison_df = pd.read_csv(comparison_mags_tablename, index_col=0)

# TODO: you can do this better...
comp_mags = {
    band: comparison_df[band][int(comparison_aperture)-1] for band in bands
}
print(comparison_df)
print("Retrieved the following comparison SDSS mags:")
print(comp_mags)
print()

target_mags = {
    band: comparison_df[band][0] for band in bands
}
print("Retrieved the following Target SDSS mags:")
print(target_mags)


##### CALCULATE K FACTORS #####

# Colour correction flux factor calculation
# I need to compute this factor:
# 10^[a({g_targ - r_targ} - {g_comp - r_comp})]
# where a is the colour term, and g-r could be u-g for each
# This factor will be close to constant across the lightcurve, so I'm
# calculating it based on the average magnitudes of the star. These
# will be sigma-clipped means!
Ku = -0.4 * a_u * ((target_mags['u'] - target_mags['g']) - (comp_mags['u'] - comp_mags['g']))
Kg = -0.4 * a_g * ((target_mags['g'] - target_mags['r']) - (comp_mags['g'] - comp_mags['r']))
Kr = -0.4 * a_r * ((target_mags['g'] - target_mags['r']) - (comp_mags['g'] - comp_mags['r']))

Ku = np.power(10.0, Ku)
Kg = np.power(10.0, Kg)
Kr = np.power(10.0, Kr)

print("\n\nFor the equation: F_targ = K * F_comp * C_targ/C_comp")
print("I calculated the following K:")
print("Ku: {:.3f}".format(Ku))
print("Kg: {:.3f}".format(Kg))
print("Kr: {:.3f}".format(Kr))


# Convert comparison magnitudes to mJy. I'll need it in a bit.
comp_flx = {key: sdss_mag2flux(comp_mags[key]) for key in comp_mags.keys()}
print("\n\nI calculated the following fluxes, in mJy: ")
print(comp_flx)

# In order to translate the GPS-stamped earth arrival times of the
# data into Barycentric MJD times, I need an astropy skycoord object.
star_loc = coord.SkyCoord(
    target_coords,
    unit=(u.hourangle, u.deg), frame='icrs'
)
print(star_loc)

# Saving data location
if not isdir(lc_dir):
    mkdir(lc_dir)
if not isdir('figs'):
    mkdir('figs')

# Loop through the logfiles, and compute the flux calibration stuff
for fname in fnames:
    # Extract the raw count data from the CCD reduction
    data = hcam.hlog.Hlog.read(fname)

    # I need to know the time for an exposure in each filter
    exptimes = {
        'r': data['1']['Exptim'].mean(),
        'g': data['2']['Exptim'].mean(),
        'u': data['3']['Exptim'].mean(),
    }
    print("\n\nThe data has exposure times of: ")
    print("r: {r:.3f}\ng: {g:.3f}\nu: {u:.3f}\n".format(**exptimes))

    print("Extracting the count fluxes...")
    target_countcurves = {
        'r': data.tseries('1', '1') / exptimes['r'],
        'g': data.tseries('2', '1') / exptimes['g'],
        'u': data.tseries('3', '1') / exptimes['u'],
    }
    comparison_countcurves = {
        'r': data.tseries('1', comparison_aperture) / exptimes['r'],
        'g': data.tseries('2', comparison_aperture) / exptimes['g'],
        'u': data.tseries('3', comparison_aperture) / exptimes['u'],
    }
    print("Done!\n")

    print("Calculating fluxes in mJy...")
    # Lets compute the fluxes of my lightcurves
    target_lightcurves = {}

    target_lightcurves['u'] = (target_countcurves['u'] / comparison_countcurves['u']) * (Ku * comp_flx['u'])
    target_lightcurves['g'] = (target_countcurves['g'] / comparison_countcurves['g']) * (Kg * comp_flx['g'])
    target_lightcurves['r'] = (target_countcurves['r'] / comparison_countcurves['r']) * (Kr * comp_flx['r'])
    print("Done!\n")

    # The flux is now sorted! But the time axis is out of whack.
    print("Sorting out the time axis (correcting to BMJD, phase folding)...")
    # Correct the MJD times recorded by the camera,
    target_lightcurves = {key: tcorrect(curve, star_loc, obsname) for key, curve in target_lightcurves.items()}

    # Phase fold. I need to know what eclipse number this one is, first.
    gband_lc = target_lightcurves['g']
    meantime = np.mean(gband_lc.t)
    E = (meantime-T0) / period
    print("The mean time of this eclipse is {:.3f}.".format(meantime))
    print("From ephemeris data, I get an eclipse Number,")
    print("  E = ({:.3f} - [T0={:.3f}]) / [P={:.5f}]".format(meantime, T0, period))
    print("  E = {:.3f}".format(E))

    E = np.rint(E)

    print("I think that the eclipse spanning from {:.3f} to {:.3f} is cycle number {}".format(
        gband_lc.t[0], gband_lc.t[-1], E)
    )

    eclTime = T0 + E*period
    print("The eclipse is then at time {:.3f}".format(eclTime))
    print("")

    # Phase fold the lightcurve
    print("Slicing out data between phase {} and {}...".format(lower_phase, upper_phase))
    for band in bands:
        # slice out the data between phase -0.5 and 0.5
        slice_time = (target_lightcurves[band].t - eclTime) / period
        slice_args = (slice_time < upper_phase)  *  (slice_time > lower_phase)

        target_lightcurves[band] = hcam.hlog.Tseries(
            slice_time[slice_args],
            target_lightcurves[band].y[slice_args],
            target_lightcurves[band].ye[slice_args],
            target_lightcurves[band].mask[slice_args]
        )

        # Bad data has error = -1
        mask = np.where(target_lightcurves[band].ye != -1)
        target_lightcurves[band] = target_lightcurves[band][mask]
    print("Calibrated file {}!".format(fname))


    # Plot the resulting lightcurve
    fig, axs = plt.subplots(3, sharex=True)

    axs[0].errorbar(
        target_lightcurves['u'].t, target_lightcurves['u'].y,
        yerr=target_lightcurves['u'].ye,
        color='blue', drawstyle='steps'
    )
    axs[0].axhline(sdss_mag2flux(target_mags['u']), color='black', label='Calculated target magnitude')
    clipped_mean_flux, _, _ = sigma_clipped_stats(target_lightcurves['u'].y)
    axs[0].axhline(clipped_mean_flux, color='magenta', label='Mean lightcurve')

    axs[1].errorbar(
        target_lightcurves['g'].t, target_lightcurves['g'].y,
        yerr=target_lightcurves['g'].ye,
        color='green', drawstyle='steps'
    )
    axs[1].axhline(sdss_mag2flux(target_mags['g']), color='black')
    clipped_mean_flux, _, _ = sigma_clipped_stats(target_lightcurves['g'].y)
    axs[1].axhline(clipped_mean_flux, color='magenta')

    axs[2].errorbar(
        target_lightcurves['r'].t, target_lightcurves['r'].y,
        yerr=target_lightcurves['r'].ye,
        color='red', drawstyle='steps'
    )
    axs[2].axhline(sdss_mag2flux(target_mags['r']), color='black')
    clipped_mean_flux, _, _ = sigma_clipped_stats(target_lightcurves['r'].y)
    axs[2].axhline(clipped_mean_flux, color='magenta')

    axs[0].set_title("Flux calibrated, phase-folded lightcurves\n{}".format(fname))
    axs[2].set_xlabel("Phase")
    axs[1].set_ylabel("SDSS Flux, mJy")

    axs[0].legend()

    date = time.Time(eclTime, format='mjd')
    date = date.strftime("%Y-%m-%d@%Hh%Mm")
    plt.savefig("figs/{}_{}.pdf".format(oname, date))
    plt.show()


    # Save the resulting lightcurves
    for key, lc in target_lightcurves.items():
        filename = "{}_{}_{}.calib".format(oname, date, key)
        filename = join(lc_dir, filename)

        with open(filename, 'w') as f:
            f.write("# Phase, Flux, Err_Flux\n")
            for t, y, ye, mask in zip(lc.t, lc.y, lc.ye, lc.mask):
                if not mask:
                    f.write("{} {} {}\n".format(t, y, ye))
