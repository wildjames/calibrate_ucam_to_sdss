import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysynphot as S
from scipy.optimize import minimize
from ucam_thruput import getref

generate_koester = True
generate_MIST = False

telescope, instrument = 'ntt', 'ucam'
stimtype = 'abmag'

variables_fname = "FOUND_VALUES.json"
if os.path.isfile(variables_fname):
    with open(variables_fname, 'r') as f:
        variables = json.load(f)
else:
    variables = {}
if not os.path.isdir("figs"):
    os.mkdir('figs')

mydir = os.path.dirname(__file__)

# These are not relevant - leave as zero for all bands
k_ext = {
    'sdss:u': 0.0,
    'sdss:g': 0.0,
    'sdss:r': 0.0,
    'sdss:i': 0.0,
    'sdss:z': 0.0,
    'ucam:ntt:u_s': 0.0,
    'ucam:ntt:g_s': 0.0,
    'ucam:ntt:r_s': 0.0,
    'ucam:ntt:i_s': 0.0,
    'ucam:ntt:z_s': 0.0,
}

def get_phoenix_mags(teff, logg, ebv):
    '''Folds a phoenix model through the SDSS lightpath, 
    and the ULTRACAM lightpath. Returns all the abmags, as surface magnitudes.'''
    row = {
        'teff': teff, 
        'logg':logg,
        'ebv': ebv,
    }

    filters = [
        'u_s', 'g_s', 'r_s', 'i_s', 'z_s', 
        # 'u', 'g', 'r', 'i', 'z'
    ]
    sdss_filters = ['u', 'g', 'r', 'i', 'z']

    sp = S.Icat('phoenix', teff, 0.0, logg)
    sp *= S.reddening.Extinction(ebv)

    # Get the *CAM telescope info
    S.setref(**getref(telescope))

    #TODO: Check this value
    airmass = 1.3
    for filt in filters:
        bp = S.ObsBandpass("{},{},{}".format(telescope,instrument,filt))
        spec = S.Observation(sp, bp, force='taper')
        mag = spec.effstim(stimtype)
        mag -= k_ext['{}:{}:{}'.format(instrument, telescope, filt)] * airmass

        row['{}:{}:{}'.format(instrument, telescope, filt)] = mag

    # Unset the *CAM thruput stuff
    S.setref(comptable=None, graphtable=None)  # leave the telescope area as it was
    # S.setref(area=None)  # reset the telescope area as well?

    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    #=-=-=-=-=   SDSS   -=-=-=-=-=#
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    airmass = 1.3
    for filt in sdss_filters:
        bp = S.ObsBandpass('sdss,{}'.format(filt))
        spec = S.Observation(sp, bp)
        mag = spec.effstim(stimtype)
        mag -= k_ext['sdss:{}'.format(filt)] * airmass

        row['sdss:{}'.format(filt)] = mag

    return row

if generate_MIST:
    # MS stars are not randomly distributed in teff/logg space. Here's a simulated isochrone.
    MIST_model = np.genfromtxt(os.path.join(mydir, "tables/MIST_ISO_log_age_8.5.csv"), delimiter=',', skip_header=True).T
    masses, teffs, loggs = MIST_model

    # Initialise the dataframe
    MIST_df = pd.DataFrame()
    # and loop through the temperatures and grvities, generating synthetic instrumental magnitudes as we go
    for logg, teff in zip(teffs, loggs):
        # Apply a temperature cut
        if teff > 4000:
            mags = get_phoenix_mags(teff, logg, 0.0)

            # I'll want these later.
            mags['u-g'] = mags['sdss:u'] - mags['sdss:g']
            mags['g-r'] = mags['sdss:g'] - mags['sdss:r']
            mags['r-i'] = mags['sdss:r'] - mags['sdss:i']
            mags['i-z'] = mags['sdss:i'] - mags['sdss:z']
            mags['g-i'] = mags['sdss:g'] - mags['sdss:i']

            for band in ['u', 'g', 'r', 'i', 'z']:
                mags["{0}-{0}_s".format(band)] = mags['{}:{}:{}_s'.format(instrument, telescope, band)] - mags['sdss:{}'.format(band)]

            MIST_df = MIST_df.append(mags, ignore_index=True)

    MIST_df.to_csv(os.path.join(mydir, "tables/MIST_mags.csv"))
else:
    MIST_df = pd.read_csv(os.path.join(mydir, "tables/MIST_mags.csv"))

################## KOESTER MODELS ##################

logg = '850'
if generate_koester:
    # Just do one gravity.
    teffs = np.arange(5000, 20001, 250)

    # SDSS filters
    sdss_filters = ['u', 'g', 'r', 'i', 'z']
    # Super SDSS filters
    cam_filters = ['u_s', 'g_s', 'r_s', 'i_s', 'z_s']

    # Lets store it all in a pandas dataframe.
    INFO = ['Teff', 'logg']
    SDSS_COLOURS = ['u-g', 'g-r', 'r-i', 'i-z', 'g-i']
    CORRECTIONS = ["{}-{}".format(a,b) for a,b in zip(sdss_filters, cam_filters)]

    COLNAMES = INFO + SDSS_COLOURS + CORRECTIONS
    koester_df = pd.DataFrame(columns=COLNAMES)

    def get_fname(teff, logg):
        'no decimal points in the logg please'
        teff = str(teff)
        logg = str(logg)
        return os.path.join(mydir, "tables/koester2/da{:>05s}_{:>03s}.dk.dat.txt".format(teff, logg))

    # I want to do each file in the koester catalogue
    for teff in teffs:
        filename = get_fname(teff, logg)
        # Unset the *CAM thruput stuff
        S.setref(comptable=None, graphtable=None)  # leave the telescope area as it was
        S.setref(area=None)  # reset the telescope area as well

        row = {
            'Teff': teff,
            'logg': float(logg)/100.0,
        }

        koester_spectrum = pd.read_csv(
            filename,
            delim_whitespace=True,
            comment='#',
            names=['WAVELENGTH (ANGSTROM)', 'FLUX (ERG/CM2/S/A)']
        )

        # drop duplicate wavelengths. WTF are they here Koester?
        koester_spectrum.drop_duplicates(subset='WAVELENGTH (ANGSTROM)', inplace=True)

        # create pysynphot spectrum
        sp = S.ArraySpectrum(
            koester_spectrum['WAVELENGTH (ANGSTROM)'],
            koester_spectrum['FLUX (ERG/CM2/S/A)'],
            fluxunits='flam',
            waveunits='Angstroms'
        )

        # Get all the SDSS magnitudes
        S.setref(comptable=None, graphtable=None)  # leave the telescope 
        simulated_mags = {}
        for f in sdss_filters:
            bp = S.ObsBandpass("{},{}".format('sdss', f))
            obs = S.Observation(sp, bp, force='taper')
            mag = obs.effstim("abmag")
            simulated_mags[f] = mag

        # Set the telescope data to the *CAM
        S.setref(**getref(telescope))

        # Magnitudes
        for f in cam_filters:
            bp = S.ObsBandpass("{},{},{}".format(telescope,instrument,f))
            obs = S.Observation(sp, bp, force='taper')
            mag = obs.effstim("abmag")
            simulated_mags[f] = mag
            row[f] = mag

        for colour in SDSS_COLOURS:
            f1, f2 = colour.split("-")
            colour_mag = simulated_mags[f1] - simulated_mags[f2]

            row[colour] = colour_mag

        # Get the actual corrections
        for colour in CORRECTIONS:
            f1, f2 = colour.split('-')
            correction_mag = simulated_mags[f1] - simulated_mags[f2]

            row[colour] = correction_mag

        # Update the table. Apply a colour cut on g-r
        koester_df = koester_df.append(row, ignore_index=True, sort=True)

    koester_df.to_csv(os.path.join(mydir, "tables/koester_magnitudes_logg_{}.csv".format(logg)))
else:
    koester_df = pd.read_csv(os.path.join(mydir, "tables/koester_magnitudes_logg_850.csv"))


################ FITTING THE TRACKS ################

targetband = input("What band are we correcting? (ugriz)  ")
diagnostic = input("What colour should we use for a diagnostic (e.g. 'u-g'): ")

def chisq(args):
    '''Uses this model to calculate chi squared:
    
    g_sdss - g_inst = g_zp + a_g(g-r)_sdss
    
    optimises for g_zp and a_g
    '''
    chisq = 0.0
    zero_point, colour_term = args

    calc_wd_correction = zero_point + (colour_term * koester_df[diagnostic])
    diff = calc_wd_correction - koester_df['{0}-{0}_s'.format(targetband)]
    chisq += (diff**2).sum()

    color = MIST_df[diagnostic]
    calc_ms_correction = zero_point + (colour_term * color)
    diff = calc_ms_correction - (MIST_df['sdss:{0}'.format(targetband)] - MIST_df['{}:{}:{}_s'.format(instrument, telescope, targetband)])
    chisq += (diff**2).sum()

    return chisq

# Run the minimiser
x0 = np.array([24.0, 0.1])
soln = minimize(
    chisq, x0
)
print("Fitting a straight line resulted in the following:")
print(soln)


fig, ax = plt.subplots(figsize=(10, 5))

MIST_df['{0}-{0}_s'.format(targetband)] = MIST_df['sdss:{0}'.format(targetband)] - MIST_df['{}:{}:{}_s'.format(instrument, telescope, targetband)]

MIST_df.plot.scatter(
    diagnostic, '{0}-{0}_s'.format(targetband),
    ax=ax,
    color='black', label="MIST MS models, 8.5Gyrs",
)
koester_df.plot.scatter(
    diagnostic, '{0}-{0}_s'.format(targetband),
    ax=ax,
    color='red', label='WD Koester Models, logg {}'.format(logg),
)

# Best fit line
zp, colterm = soln['x']

xr = np.linspace(MIST_df[diagnostic].min(), MIST_df[diagnostic].max(), 100)
yr = zp + colterm * xr

ax.plot(xr, yr, color='blue', label='Fitted line, colour term: {:.3f}'.format(colterm))

ax.set_ylabel("SDSS ${}$ - {}/{} ${}_s$".format(targetband, telescope.upper(), instrument.upper(), targetband))
ax.set_xlabel("SDSS ${}$".format(diagnostic))


# SA 114 SDSS colour
if diagnostic == 'g-r':
    ax.axvline(1.120, color='magenta', label='SA 114 548')
else:
    ax.axvline(3.146, color='magenta', label='SA 114 548')


ax.legend()
plt.tight_layout()
plt.savefig("figs/{}band.pdf".format(targetband))
plt.show()

# Save the new value to file
variables["a{}".format(targetband)] = colterm
with open(variables_fname, 'w') as f:
    f.write(json.dumps(variables))
