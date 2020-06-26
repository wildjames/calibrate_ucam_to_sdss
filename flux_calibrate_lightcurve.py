from pprint import pprint
from time import sleep

import hipercam as hcam
import matplotlib.pyplot as plt
import numpy as np

from calphot.constructReference import get_instrumental_mags

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
########## USER DEFINED INPUTS ##########
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
# Where are the target data stored? 
# Aperture 1 is assumed to be the target star!
fname = 'data/run015.log'
target_coords = "20 29 17.13 -43 40 19.8"

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

fig, axs = plt.subplots(3, sharex=True)

axs[0].errorbar(target_lightcurves['u'].t, target_lightcurves['u'].y, yerr=target_lightcurves['u'].ye, color='blue', drawstyle='steps')
axs[1].errorbar(target_lightcurves['g'].t, target_lightcurves['g'].y, yerr=target_lightcurves['g'].ye, color='green', drawstyle='steps')
axs[2].errorbar(target_lightcurves['r'].t, target_lightcurves['r'].y, yerr=target_lightcurves['r'].ye, color='red', drawstyle='steps')

axs[2].set_ylabel("Time, MJD")
axs[1].set_ylabel("SDSS Flux, mJy")

plt.show()
