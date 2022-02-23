import json
from os.path import abspath, join, split
from time import sleep

import hipercam as hcam
import pandas as pd

from constructReference import get_instrumental_mags

# Gather some data files
mydir = split(abspath(__file__))[0]
# This is the machine-readable version of the smith table. You shouldn't need to fiddle this!
smith_table_fname = join(mydir, "tables/tab08.dat.txt")
values_fname = "FOUND_VALUES.json"
print("Using the smith (2002) data found in:")
print(smith_table_fname)
print("and the prior variables stored at: ")
print(abspath(values_fname))
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

    print("\n\nGrabbing the magnitudes of {} from the smith table at {}".format(starname, tablename))
    # We gotta build it by hand
    gathered_mags = {}
    gathered_mags['r'] = std_data["r'"].values[0]
    gathered_mags['g'] = gathered_mags["r"] + std_data["g'-r'"].values[0]
    gathered_mags['u'] = gathered_mags["g"] + std_data["u'-g'"].values[0]
    gathered_mags['i'] = gathered_mags["r"] - std_data["r'-i'"].values[0]
    gathered_mags['z'] = gathered_mags["i"] - std_data["i'-z'"].values[0]

    coords = "{} {}".format(std_data['RA (J2000.0)'].values[0], std_data['DEC (J2000.0)'].values[0])
    return gathered_mags, coords


def gather_standard_mags(fname, coords, obsname, k_ext, bands):
    '''return the instrumental, above atmosphere, SDSS magnitudes.
    the passed list of bands tells me what order the CCDs are in.'''
    data = hcam.hlog.Hlog.rascii(fname)
    inst_mags = get_instrumental_mags(data, coords, obsname, k_ext)
    inst_mags = {band: inst_mags[str(i+1)][0] for i, band in enumerate(bands)}
    return inst_mags


if __name__ in "__main__":
    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#
    #=--=#=--= USER DEFINED STUFF =--=#=--=#
    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#
    desc = '''
    The standard star reduction MUST be done with large apertures to
    capture all the light from the star.

    The Smith tables contain good, reliable magnitudes for the standards.
    These can be retrieved by entering its name below, but check the format is correct first!


    NOTE!!! This script assumes that the colour terms are all a
    function of g-r, EXCEPT the u band, which is assumed to be a
    function of u-g. IF YOU CHANGE THAT, CHANGE THIS SCRIPT!!!
    '''
    import argparse

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('std_reduction')
    parser.add_argument('std_name')

    # standard_star_reduction = "Quality_Reductions/2019_09_27-run007.log"
    # std_name = 'SA 114 548'

    args = parser.parse_args()
    standard_star_reduction = args.std_reduction
    std_name = args.std_name

    # Observed bands in the above reduction.
    bands = input("Please enter the filters, in order, separated by spaces: ").split(" ")
    print(" -->> I have the bands: {}".format(bands))
    # Observing conditions
    obsname = input("Please enter the observation site: ")
    # Atmospheric extinction, in CCD order (typically r,g,u)
    k_ext = [0.1129, 0.2020, 0.4868]
    k_ext = []
    for band in bands:
        k_ext.append(float(input("Please enter the extinction coefficient for band {}: ".format(band))))

    print(" -->> I have the extinction coefficients:")
    for band, k in zip(bands, k_ext):
        print("{:>5s}: {}".format(band, k))

    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#
    #=--=#=--=#=--=#=--=#=--=#=--=#=--=#=--=#


    # Gather data on the standard star
    smith_mags, smith_coords = get_smith_info(std_name, smith_table_fname)
    print("{} magnitudes:".format(std_name))
    print('\n'.join(["{}: {:.3f}".format(key, val) for key, val in smith_mags.items()]))

    print("I will gather the instrumental magnitudes observed in {}\n\n".format(standard_star_reduction))
    sleep(2)
    inst_mags = gather_standard_mags(standard_star_reduction, smith_coords, obsname, k_ext, bands)

    print("\n\nI have previously found these values:")
    print('\n'.join(["{}: {}".format(key, val) for key, val in variables.items()]))
    print('\n\n\n')

    print("\n\nTo calculate the zero points, I use equations of the form:")
    print("g_zp = g_inst,noatmos - g_sdss - (a_g * (g-r)_sdss)")
    zero_points = {}
    for band in bands:
        colour_term = variables["a{}".format(band)]

        print("\n{}_zp".format(band))
        if 'u' in band:
            print("  using u-g")
            colour_term *= smith_mags['u'] - smith_mags['g']
        elif 'i' in band:
            print("  using g-i")
            colour_term *= smith_mags['g'] - smith_mags['i']
        else:
            print("  using g-r")
            colour_term *= smith_mags['g'] - smith_mags['r']

        print("Colour term: {:.3f}".format(colour_term))
        zero_points[band] = smith_mags[band] - inst_mags[band] - colour_term
        variables["{}_zp".format(band)] = zero_points[band]
        print("{}_zp: {:.3f}".format(band, zero_points[band]))

    print("\n\nI calculated the zero points, accounting for airmass in the target reduction and the bandpass differences between UCAM and SDSS, as:")
    print('\n'.join(["{}: {:.4f}".format(key, val) for key, val in zero_points.items()]))

    print("Dumping values into {}".format(values_fname))
    with open(values_fname, 'w') as f:
        f.write(json.dumps(variables))
