from os.path import join, isdir
from os import mkdir
from pprint import pprint
from time import sleep

import hipercam as hcam
import matplotlib.pyplot as plt
import numpy as np
from astropy import time, coordinates as coord
from astropy import units as u


from calphot.constructReference import get_instrumental_mags
from calphot.getEclipseTimes import tcorrect

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
########## USER DEFINED INPUTS ##########
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
# Where are the target data stored? 
# Aperture 1 is assumed to be the target star!
#fname = 'data/run015.log'
fname = 'Quality_Reductions/2019_09_27-run015.log'
oname = 'ASASSN-17jf'
lc_dir = 'LIGHTCURVES'

target_coords = "20 29 17.13 -43 40 19.8"
T0 = 58754.12003
period = 0.056789
lower_phase, upper_phase = -0.5, 0.5

# Comparison star SDSS magnitudes 
#       => Target will be in SDSS mags!
comp_mags = {
    'u': 19.728,
    'g': 16.874,
    'r': 15.834,
}

# Atmospheric extinction coefficients. Calculated empirically for UCAM, NTT, super filters.
# Must be in CCD order of the logfile. This is typically r, g, u
k_ext = [0.1129, 0.2020, 0.4868]
# Observatory name
obsname = 'lasilla'

# Colour term in each band, found by 
# generate_colourtracks.py (READ THAT BEFORE YOU USE IT!!)
##### --->> a_u uses u-g, g and r are g-r. <<--- #####
a_u = -0.037
a_g = -0.024
a_r = -0.032

# Zero points in each CCD for UCAM on the NTT.
# Calculated from SA 114 548 on 27 Sept. 2019
zp_u = 24.817
zp_g = 26.218
zp_r = 25.785

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#


data = hcam.hlog.Hlog.read(fname)

star_loc = coord.SkyCoord(
    target_coords,
    unit=(u.hourangle, u.deg), frame='icrs'
)


def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux


# Convert to mJy
comp_flx = {key: sdss_mag2flux(comp_mags[key]) for key in comp_mags.keys()}

print("\n\nUsing the following SDSS magnitudes: ")
pprint(comp_mags)
print("I calculated the following fluxes, in mJy: ")
pprint(comp_flx)


# I need to know the time for an exposure in each filter
exptimes = {
    'r': data['1']['Exptim'].mean(),
    'g': data['2']['Exptim'].mean(),
    'u': data['3']['Exptim'].mean(),
}
print("\nThe data has exposure times of: ")
print("r: {r:.3f}\ng: {g:.3f}\nu: {u:.3f}\n".format(**exptimes))
sleep(3)

# Extract the raw count data from the CCD reduction
target_countcurves = {
    'r': data.tseries('1', '1') / exptimes['r'],
    'g': data.tseries('2', '1') / exptimes['g'],
    'u': data.tseries('3', '1') / exptimes['u'],
}
comparison_countcurves = {
    'r': data.tseries('1', '2') / exptimes['r'],
    'g': data.tseries('2', '2') / exptimes['g'],
    'u': data.tseries('3', '2') / exptimes['u'],
}

# Lets calcualate the instrumental magnitudes here. This subtracts the atmosphere as well.
print("Fetching instrumental magnitudes, corrected for atmospheric extinction....")
sleep(2)
target_instmags = get_instrumental_mags(data, target_coords, obsname, k_ext)
print("\n\n\n\n\nDone fetching the atmosphere-subtracted instrumental magnitudes!")

# massage the above into the correct form for this script
target_instmags = {
    'r': target_instmags['1'][0],
    'g': target_instmags['2'][0],
    'u': target_instmags['3'][0],
}

print("\nI got an average instrumental magnitude of: ")
pprint(target_instmags)
print("\nAdding the zero points of:")
print("r: {:.3f}\ng: {:.3f}\nu: {:.3f}".format(zp_r, zp_g, zp_u))

target_instmags['r'] += zp_r
target_instmags['g'] += zp_g
target_instmags['u'] += zp_u

print("\ngives uncorrected instrumental magnitudes of: ")
pprint(target_instmags)

# Initial SDSS mags are just set to the 
# instrumental, zero point-added, atmos-subtracted magnitudes to start
target_sdssmags = {
    key: val for key, val in target_instmags.items()
}

# I need to compute this factor:
# 10^[0.4*a({g_targ - r_targ} - {g_comp - r_comp})]
# where a is the colour term, and g-r could be u-g for each
# This factor will be close to constant across the lightcurve, so I'm 
# calculating it based on the average magnitudes of the star. These 
# will be sigma-clipped means! But it means I need to find the 
# SDSS colours of the target before I can continue. 
#
# Lets do that.

# Tracking variables:
du = 99
dg = 99
dr = 99

iteration = 0

while (du + dg + dr) > 0.001:
    iteration += 1

    u_sdss_new = target_instmags['u'] + a_u*(comp_mags['u'] - comp_mags['g'])
    g_sdss_new = target_instmags['g'] + a_g*(comp_mags['g'] - comp_mags['r'])
    r_sdss_new = target_instmags['r'] + a_r*(comp_mags['g'] - comp_mags['r'])

    du = abs(target_sdssmags['u'] - u_sdss_new)
    dg = abs(target_sdssmags['g'] - g_sdss_new)
    dr = abs(target_sdssmags['r'] - r_sdss_new)

    target_sdssmags['u'] = u_sdss_new
    target_sdssmags['g'] = g_sdss_new
    target_sdssmags['r'] = r_sdss_new

    print("Iteration {:>03d} | du: {:>06.3f} | dg: {:>06.3f} | dr: {:>06.3f} |".format(iteration, du, dg, dr), end='\n')

print("\n\n\nConverged on the following target SDSS magnitudes:")
print("r: {r:.3f}\ng: {g:.3f}\nu: {u:.3f}\n".format(**target_sdssmags))

colour_terms = {}
colour_terms['u'] = a_u*(comp_mags['u'] - comp_mags['g'])
colour_terms['g'] = a_g*(comp_mags['g'] - comp_mags['r'])
colour_terms['r'] = a_r*(comp_mags['g'] - comp_mags['r'])
print("The colour corrections applied to this star are:")
print("r: {r:.3f}\ng: {g:.3f}\nu: {u:.3f}\n".format(**colour_terms))

# What do these correspond to?
k_u = 0.4 * a_u * (target_sdssmags['u'] - target_sdssmags['g'] - comp_mags['u'] + comp_mags['g'])
k_g = 0.4 * a_g * (target_sdssmags['g'] - target_sdssmags['r'] - comp_mags['g'] + comp_mags['r'])
k_r = 0.4 * a_r * (target_sdssmags['g'] - target_sdssmags['r'] - comp_mags['g'] + comp_mags['r'])

k_u = np.power(10,k_u)
k_g = np.power(10,k_g)
k_r = np.power(10,k_r)

print("The following are the K factors in the equation: F_targ = K * F_comp * C_targ/C_comp")
print("k_r: {:.3f}\nk_g: {:.3f}\nk_u: {:.3f}\n".format(k_r, k_g, k_u))


target_lightcurves = {}

target_lightcurves['u'] = (target_countcurves['u'] / comparison_countcurves['u']) * (k_u * comp_flx['u'])
target_lightcurves['g'] = (target_countcurves['g'] / comparison_countcurves['g']) * (k_g * comp_flx['g'])
target_lightcurves['r'] = (target_countcurves['r'] / comparison_countcurves['r']) * (k_r * comp_flx['r'])


# Correct the MJD times recorded by the camera, 
target_lightcurves = {key: tcorrect(curve, star_loc, obsname) for key, curve in target_lightcurves.items()}


gband_lc = target_lightcurves['g']
meantime = np.mean(gband_lc.t)
E = (meantime-T0) / period
print("  The mean time of this eclipse is {:.3f}.".format(meantime))
print("  From ephemeris data, I get an eclipse Number,")
print("    E = ({:.3f} - [T0={:.3f}]) / [P={:.5f}]".format(meantime, T0, period))
print("    E = {:.3f}".format(E))

E = np.rint(E)
# The above can be off, if the eclipse isnt the minimum. in/decriment until it's within bounds
while T0 + E*period < gband_lc.t[0]:
    print("    !!! Eclipse time not within these data! Incrimenting E...")
    E += 1
while T0 + E*period > gband_lc.t[-1]:
    print("    !!! Eclipse time not within these data! Decrimenting E...")
    E -= 1

print("  I think that the eclipse spanning from {:.3f} to {:.3f} is cycle number {}".format(
    gband_lc.t[0], gband_lc.t[-1], E)
)

eclTime = T0 + E*period
print("  The eclipse is then at time {:.3f}".format(eclTime))
print("")

# Phase fold the lightcurve
print("Slicing out data between phase {} and {}".format(lower_phase, upper_phase))
for key in ['u','g','r']:
    # slice out the data between phase -0.5 and 0.5
    slice_time = (target_lightcurves[key].t - eclTime) / period
    slice_args = (slice_time < upper_phase)  *  (slice_time > lower_phase)

    target_lightcurves[key] = hcam.hlog.Tseries(
        slice_time[slice_args],
        target_lightcurves[key].y[slice_args],
        target_lightcurves[key].ye[slice_args],
        target_lightcurves[key].mask[slice_args]
    )
    
    # Bad data has error = -1
    mask = np.where(target_lightcurves[key].ye != -1)
    target_lightcurves[key] = target_lightcurves[key][mask]
print("Done!")





fig, axs = plt.subplots(3, sharex=True)

axs[0].errorbar(
    target_lightcurves['u'].t, target_lightcurves['u'].y, 
    yerr=target_lightcurves['u'].ye, 
    color='blue', drawstyle='steps'
)
axs[1].errorbar(
    target_lightcurves['g'].t, target_lightcurves['g'].y, 
    yerr=target_lightcurves['g'].ye, 
    color='green', drawstyle='steps'
)
axs[2].errorbar(
    target_lightcurves['r'].t, target_lightcurves['r'].y, 
    yerr=target_lightcurves['r'].ye, 
    color='red', drawstyle='steps'
)

axs[0].set_title("Flux calibrated, phase-folded lightcurves")
axs[2].set_xlabel("Phase")
axs[1].set_ylabel("SDSS Flux, mJy")

plt.show()



# Saving data
if not isdir(lc_dir):
    mkdir(lc_dir)

for key, lc in target_lightcurves.items():
    date = time.Time(eclTime, format='mjd')
    date = date.strftime("%Y-%m-%d@%Hh%Mm")

    filename = oname
    filename = "{}_{}_{}.calib".format(filename, date, key)
    filename = join(lc_dir, filename)

    with open(filename, 'w') as f:
        f.write("# Phase, Flux, Err_Flux\n")
        for t, y, ye, mask in zip(lc.t, lc.y, lc.ye, lc.mask):
            if not mask:
                f.write("{} {} {}\n".format(t, y, ye))
