import hipercam as hcam
import numpy as np
from pprint import pprint
from calphot.constructReference import get_instrumental_mags

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
########## USER DEFINED INPUTS ##########
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
# Where are the target data stored? 
# Aperture 1 is assumed to be the target star!
fname = 'data/run015.log'
target_coords = "20 29 17.13 -43 40 19.8"

# Colour term in each band, found by 
# generate_colourtracks.py (READ THAT BEFORE YOU USE IT!!)
##### --->> a_u uses u-g, g and r are g-r. <<--- #####
a_u = -0.037
a_g = -0.024
a_r = -0.032

# Zero points in each CCD for UCAM on the NTT.
# Calculated from SA 114 548 on 27 Sept. 2019
u_zp = 24.817
g_zp = 26.218
r_zp = 25.785

# Comparison star SDSS magnitudes 
#       => Target will be in SDSS mags!
comp_mags = {
    'u': 19.728,
    'g': 16.874,
    'r': 15.834,
}

# Atmospheric extinction coefficients. Calculated empirically for UCAM, NTT, super filters.
k_ext = [0.1129, 0.2020, 0.4868]
# Observatory name
obsname = 'lasilla'

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#


def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux


comp_flx = {key: sdss_mag2flux(comp_mags[key]) for key in comp_mags.keys()}

print("Using the following SDSS magnitudes: ")
pprint(comp_mags)
print("I calculated the following fluxes, in mJy: ")
pprint(comp_flx)

data = hcam.hlog.Hlog.read(fname)

# Extract the raw count data from the CCD reduction
target_countcurves = {
    'r': data.tseries('1', '1'),
    'g': data.tseries('2', '1'),
    'u': data.tseries('3', '1'),
}
comparison_lightcurves = {
    'r': data.tseries('1', '2'),
    'g': data.tseries('2', '2'),
    'u': data.tseries('3', '2'),
}

# Lets calcualate the instrumental magnitudes here. This subtracts the atmosphere as well.
target_instmags = get_instrumental_mags(data, target_coords, obsname, k_ext)
# massage the above into the correct form for this script
target_instmags = {
    'r': target_instmags['1'][0],
    'g': target_instmags['2'][0],
    'u': target_instmags['3'][0],
}

target_sdssmags = {
    key: val for key, val in target_instmags.iter()
}

print("I got an average instrumental magnitude of: ")
pprint(target_instmags)
print("Adding the zero points of:")
print("r: {:.3f}\ng: {:.3f}\nu: {:.3f}".format(zp_r, zp_g, zp_u))

target_instmags['r'] += zp_r
target_instmags['g'] += zp_g
target_instmags['u'] += zp_u

print("gives uncorrected instrumental magnitudes of: ")
pprint(target_instmags)