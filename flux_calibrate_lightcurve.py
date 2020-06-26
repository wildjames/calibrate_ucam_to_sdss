import hipercam as hcam
import numpy as np
from pprint import pprint


#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
########## USER DEFINED INPUTS ##########
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
# Where are the target data stored? 
# Aperture 1 is assumed to be the target star!
fname = 'data/run015.log'

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

print(dir(target_countcurves['u'].to_mag()))
target_instmags = {
    'u': target_countcurves['u'].to_mag()
}


