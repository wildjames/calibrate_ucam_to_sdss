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

from calphot.constructReference import get_instrumental_mags
from calphot.getEclipseTimes import tcorrect

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
########## USER DEFINED INPUTS ##########
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

# Where are the target data stored? 
# Aperture 1 is assumed to be the target star!
fname = 'Quality_Reductions/2019_09_27-run015.log'
oname = 'ASASSN-17jf'
lc_dir = 'LIGHTCURVES'

target_coords = "20 29 17.13 -43 40 19.8"
T0 = 58754.12003
period = 0.056789
lower_phase, upper_phase = -0.5, 0.5

# Atmospheric extinction coefficients. Calculated empirically for UCAM, NTT, super filters.
# Must be in CCD order of the logfile. This is typically r, g, u
k_ext = [0.1129, 0.2020, 0.4868]
# Observatory name
obsname = 'lasilla'

# Comparison stars data
comparison_aperture = 5
comparison_mags_tablename = 'tables/comparison_star_sdss_mags.csv'

# Previously calculated values
variables_fname = "FOUND_VALUES.json"
with open(variables_fname, 'r') as f:
    variables = json.load(f)

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=##=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=##=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=##=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

comparison_df = pd.from_csv(comparison_mags_tablename)

comp_mags = dict(comparison_df[comparison_df['aperture'] == comparison_aperture])
print("Retrieved the following comparison mags:")
print(comp_mags)

# Colour terms. a_u is u-g, g and r are g-r.
a_u = variables['au']
a_g = variables['ag']
a_r = variables['ar']

# Zero points in each CCD for UCAM on the NTT.
u_zp = variables['u_zp']
g_zp = variables['g_zp']
r_zp = variables['r_zp']



def sdss_mag2flux(mag):
    '''Takes an SDSS magnitude, returns the corresponding flux in [mJy]'''
    alpha = 3631e3

    flux = 10**(-mag/2.5)
    flux*= alpha

    return flux




data = hcam.hlog.Hlog.read(fname)

star_loc = coord.SkyCoord(
    target_coords,
    unit=(u.hourangle, u.deg), frame='icrs'
)

# Convert to mJy
comp_flx = {key: sdss_mag2flux(comp_mags[key]) for key in comp_mags.keys()}

print("\n\nUsing the following comparison SDSS magnitudes (AP{}): ".format(comparison_aperture))
pprint(comp_mags)
print("I calculated the following fluxes, in mJy: ")
pprint(comp_flx)
