import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S
from ucam_thruput import getref
import pandas as pd
import os
from pathlib import Path
from scipy.optimize import minimize


telescope, instrument = 'ntt', 'ucam'
stimtype = 'abmag'

filters = [
    'u_s', 'g_s', 'r_s', 'i_s', 'z_s', 
    # 'u', 'g', 'r', 'i', 'z'
]
sdss_filters = ['u', 'g', 'r', 'i', 'z']
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

targetband = 'u'
# diagnostic = 'u-g'
diagnostic = 'g-r'

def get_all_mags(teff, logg, ebv):
    '''Folds a phoenix model through the SDSS lightpath, and the ULTRACAM lightpath. Returns all the abmags, as surface magnitudes.'''
    row = {
        'teff': teff, 
        'logg':logg,
        'ebv': ebv,
    }

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


MIST_model = np.genfromtxt("tables/MIST_ISO_log_age_8.5.csv", delimiter=',', skip_header=True).T
masses, teffs, loggs = MIST_model

MIST_df = pd.DataFrame()
for logg, teff in zip(teffs, loggs):
    if teff > 4000:
        mags = get_all_mags(teff, logg, 0.0)

        mags['g-r'] = mags['sdss:g'] - mags['sdss:r']
        mags['u-g'] = mags['sdss:u'] - mags['sdss:g']

        mags["{0}_s-{0}".format(targetband)] = mags['ucam:ntt:{}_s'.format(targetband)] - mags['sdss:{}'.format(targetband)]

        MIST_df = MIST_df.append(mags, ignore_index=True)

MIST_df.to_csv("MIST_mags.csv")

################## KOESTER MODELS ##################

logg = '850'
teffs = [
    5000, 5250, 5500, 5750, 
    6000, 6250, 6500, 6750, 
    7000, 7250, 7500, 7750, 
    8000, 8250, 8500, 8750, 
    9000, 9250, 9500, 9750, 
    10000, 10250, 10500, 10750, 
    11000, 11250, 11500, 11750, 
    12000, 12250, 12500, 12750, 
    13000, 13250, 13500, 13750, 
    13000, 13250, 13500, 13750, 
    13000, 13250, 13500, 13750, 
    14000, 14250, 14500, 14750, 
    15000, 15250, 15500, 15750, 
    16000, 16250, 16500, 16750, 
    17000, 17250, 17500, 17750, 
    18000, 18250, 18500, 18750, 
    19000, 19250, 19500, 19750, 
    20000
]

# SDSS filters
sdss_filters = ['u', 'g', 'r', 'i', 'z']
# Super SDSS filters
cam_filters = ['u_s', 'g_s', 'r_s', 'i_s', 'z_s']

# Lets store it all in a pandas dataframe.
INFO = ['Teff', 'logg']
CAM_COLOURS = ['u_s-g_s', 'g_s-r_s', 'r_s-i_s', 'i_s-z_s']
SDSS_COLOURS = ['u-g', 'g-r', 'r-i', 'i-z']
CORRECTIONS = ["{}-{}".format(a,b) for a,b in zip(cam_filters, sdss_filters)]

COLNAMES = INFO + SDSS_COLOURS + CAM_COLOURS + CORRECTIONS
koester_df = pd.DataFrame(columns=COLNAMES)

def get_fname(teff, logg):
    'no decimal points in the logg please'
    teff = str(teff)
    logg = str(logg)
    return "koester2/da{:>05s}_{:>03s}.dk.dat.txt".format(teff, logg)

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

    # Get the colour corrections for the super filters
    S.setref(**getref(telescope))

    # Magnitudes
    for f in cam_filters:
        bp = S.ObsBandpass("{},{},{}".format(telescope,instrument,f))
        obs = S.Observation(sp, bp, force='taper')
        mag = obs.effstim("abmag")
        simulated_mags[f] = mag

    # Colours
    for colour in CAM_COLOURS:
        f1, f2 = colour.split('-')
        correction_mag = simulated_mags[f1] - simulated_mags[f2]

        row[colour] = correction_mag

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

koester_df.to_csv("koester_magnitudes_logg_{}.csv".format(logg))


def chisq(args):
    chisq = 0.0
    zero_point, colour_term = args

    calc_wd_correction = zero_point + colour_term * koester_df[diagnostic]

    diff = calc_wd_correction - koester_df['{0}_s-{0}'.format(targetband)]
    chisq += (diff**2).sum()


    if diagnostic == 'u-g':
        color = MIST_df['sdss:u']-MIST_df['sdss:g']
    else:
        color = MIST_df['sdss:g']-MIST_df['sdss:r']
    calc_ms_correction = zero_point + (colour_term * color)

    diff = calc_ms_correction - (MIST_df['ucam:ntt:{0}_s'.format(targetband)] - MIST_df['sdss:{0}'.format(targetband)])
    chisq += (diff**2).sum()

    return chisq


x0 = np.array([24.0, 0.1])
soln = minimize(
    chisq, x0
)
print(soln)


fig, ax = plt.subplots(figsize=(10, 5))

MIST_df.plot.scatter(
    diagnostic, '{0}_s-{0}'.format(targetband),
    ax=ax,
    color='black', label="MIST MS models, 8.5Gyrs",
)
koester_df.plot.scatter(
    diagnostic, '{0}_s-{0}'.format(targetband),
    ax=ax,
    color='red', label='WD Koester Models, logg {}'.format(logg),
)


zp, colterm = soln['x']

xr = np.linspace(MIST_df[diagnostic].min(), MIST_df[diagnostic].max(), 100)
yr = zp + colterm * xr

ax.plot(xr, yr, color='blue', label='Fitted line, colour term: {:.3f}'.format(colterm))

# ax.set_xlabel("UCAM/NTT gs - rs")
ax.set_ylabel("UCAM/NTT ${0}_s$ - SDSS ${0}$".format(targetband))
ax.set_xlabel("SDSS ${}$".format(diagnostic))


# Data we care about. This should be in the same colour as the tables!
# # UCAM/NTT instrumental g_s-r_s
# # ASASSN-17jf 
# ax.axvline(-0.701918, color='orange', label='ASASSN-17jf')

# # Comparison stars
# ax.axvline( 0.598594, color='blue')
# ax.axvline(-0.042797, color='blue')
# ax.axvline( 0.056079, color='blue')
# ax.axvline( 0.702917, color='blue')
# ax.axvline( 0.233909, color='blue')
# ax.axvline(-0.094819, color='blue', label='Comparison Stars')

# SA 114 SDSS colour
if diagnostic == 'g-r':
    ax.axvline(1.120, color='magenta', label='SA 114 548')
else:
    ax.axvline(3.146, color='magenta', label='SA 114 548')


ax.legend()
plt.tight_layout()
plt.savefig("figs/{}band.pdf".format(targetband))
plt.show()

