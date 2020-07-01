import json

# import hipercam as hcam
import numpy as np
import pandas as pd
from os.path import abspath, split, join

# from calphot.constructReference import get_instrumental_mags

mydir = split(abspath(__file__))[0]
# This is the machine-readable version of the smith table. You shouldn't need to fiddle this!
smith_table_fname = join(mydir, "tables/tab08.dat.txt")
values_fname = join(mydir, "FOUND_VALUES.json")
print("Using the smith (2002) data found in:")
print(smith_table_fname)
print("and the prior variables stored at: ")
print(values_fname)
# First, remember what we've found so far
with open(values_fname, 'r') as f:
    variables = json.load(f)

def get_smith_info(starname, tablename):
    '''Gather the information I care about from the smith data table. 
    The smith star name can be finnicky, if this doesn't work, check how the star
    is formatted in the table.
    
    Returns:
    --------
    gathered_mags, dict:
        The magnitudes of the star, keyed by filter
    coords, str
        The coordinates of the star, as a string formatted in a way readable by astropy.
    '''
    # Load up the smith data, and grab our star
    smith_table = pd.read_fwf(tablename)
    std_data = smith_table.loc[smith_table['StarName']==starname]

    print("Grabbing the magnitudes of {} from the smith table at {}".format(starname, tablename))
    # We gotta build it by hand
    gathered_mags = {}
    gathered_mags['r'] = std_data["r'"].values[0]
    gathered_mags['g'] = gathered_mags["r"] + std_data["g'-r'"].values[0]
    gathered_mags['u'] = gathered_mags["g"] + std_data["u'-g'"].values[0]
    gathered_mags['i'] = gathered_mags["r"] - std_data["r'-i'"].values[0]
    gathered_mags['z'] = gathered_mags["i"] - std_data["i'-z'"].values[0]

    coords = "{} {}".format(std_data['RA (J2000.0)'].values[0], std_data['RA (J2000.0)'].values[0])
    return gathered_mags, coords


def gather_standard_mags(fname, coords, obsname, k_ext, bands):
    '''return the instrumental, above atmosphere, SDSS magnitudes. 
    the passed list of bands tells me what order the CCDs are in.'''
    inst_mags = get_instrumental_mags(fname, coords, obsname, k_ext)
    inst_mags = {band: inst_mags[str(i+1)] for i, band in enumerate(bands)}
    #     'r': inst_mags['1'][0],
    #     'g': inst_mags['2'][0],
    #     'u': inst_mags['3'][0],
    # }
    return inst_mags


if __name__ in "__main__":
    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#
    #=--=#=--= USER DEFINED STUFF =--=#=--=#
    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#
    '''
    The standard star reduction MUST be done with large apertures to 
    capture all the light from the star.

    The Smith tables contain good, reliable magnitudes for the standards. 
    These can be retrieved by entering its name below, but check the format is correct first!


    NOTE!!! This script assumes that the colour terms are all a 
    function of g-r, EXCEPT the u band, which is assumed to be a 
    function of u-g. IF YOU CHANGE THAT, CHANGE THIS SCRIPT!!!
    '''

    standard_star_reduction = "run007.log"
    std_name = 'SA 114 548'

    # Observed bands in the above reduction.
    bands = ['r', 'g', 'u']
    # Observing conditions
    obsname = 'lasilla'
    # Atmospheric extinction, in CCD order (typically r,g,u)
    k_ext = [0.1129, 0.2020, 0.4868]

    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#
    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#


    print("I have previously found these values:")
    print('\n'.join(["{}: {}".format(key, val) for key, val in variables.items()]))
    print('\n\n\n')

    # Gather data on the standard star
    smith_mags, smith_coords = get_smith_info(std_name, smith_table_fname)
    print("{} magnitudes:".format(std_name))
    print('\n'.join(["{}: {}".format(key, val) for key, val in smith_mags.items()]))

    inst_mags = gather_standard_mags(standard_star_reduction, smith_coords, obsname, k_ext, bands)

    zero_points = {}
    for band in bands:
        colour_term = variables["a{}".format(band)]

        if 'u' in band:
            colour_term *= smith_mags['u'] - smith_mags['g']
        else:
            colour_term *= smith_mags['g'] - smith_mags['r']

        zero_points[band] = inst_mags[band] - smith_mags[band]
        variables["{}_zp".format(band)] = zero_points[band]

    print("I calculated the zero points, accounting for airmass in the target reduction and the bandpass differences between UCAM and SDSS, as:")
    print('\n'.join(["{}: {}".format(key, val) for key, val in zero_points.items()]))

    with open(values_fname, 'w') as f:
        f.write(json.dumps(variables))
