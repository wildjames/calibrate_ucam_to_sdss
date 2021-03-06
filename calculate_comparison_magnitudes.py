import argparse
import json
from os.path import abspath, join, split

import hipercam as hcam
import pandas as pd

from calphot.constructReference import get_instrumental_mags

desc = '''
Using the standard - like reduction, calculate the comparison star magnitudes of a logfile.
g_sdss = g_inst,noatmos - g_zp - a_g(g-r) 
'''


parser = argparse.ArgumentParser(description=desc, prefix_chars='@')
parser.add_argument("reduction")
parser.add_argument("RA")
parser.add_argument("DEC")
parser.add_argument('@@oname', default="comparison_star_sdss_mags", required=False)

args = parser.parse_args()

print(args)

fname = args.reduction
target_coords = "{} {}".format(args.RA, args.DEC)
oname = args.oname

obsname = 'lasilla'
k_ext = [0.1129, 0.2020, 0.4868]
bands = ['r', 'g', 'u']

obsname = input("Please enter the observation site name: ")

# Gather some data files
# mydir = split(abspath(__file__))[0]
values_fname = abspath("FOUND_VALUES.json")

print("and the prior variables stored at: ")
print(values_fname)

# First, remember what we've found so far
with open(values_fname, 'r') as f:
    variables = json.load(f)


data = hcam.hlog.Hlog.read(fname)
# Gets a dict of lists, of all the aperture mags. The first one will be the target, so can be discarded.
inst_mags = get_instrumental_mags(data, target_coords, obsname, k_ext)
inst_mags = {bands[int(key)-1]: mags for key, mags in inst_mags.items()}

print("I got the following comparison instrumental magnitudes:")
for key, mags in inst_mags.items():
    print("Band: {}".format(key))
    for mag in mags:
        print("  {:.3f}".format(mag))
    print()


# Compute the SDSS magnitude of each aperture. You need to do this iteratively, and 
# simultaneously in all bands, since the colour term used to calculate the SDSS magnitude 
# from instrumental DEPENDS on the SDSS magnitude!
sdss_mags = {band: [] for band in bands}

# Colour terms. a_u is u-g, g and r are g-r.
a_u = variables['au']
a_g = variables['ag']
a_i = variables['ar']

# Zero points in each CCD for UCAM on the NTT.
u_zp = variables['u_zp']
g_zp = variables['g_zp']
i_zp = variables['r_zp']

for ap_index in range(len(mags)):
    u_inst = inst_mags['u'][ap_index]
    g_inst = inst_mags['g'][ap_index]
    i_inst = inst_mags['r'][ap_index]
    print("Got the following instrumental magnitudes for aperture {}:".format(ap_index+1))
    print("r: {:.3f}".format(r_inst))
    print("g: {:.3f}".format(g_inst))
    print("u: {:.3f}".format(u_inst))

    # First iteration uses no colour term, so set the sdss 
    # mags as equal to the instrumental ones.
    u_sdss = u_inst + u_zp
    g_sdss = g_inst + g_zp
    r_sdss = r_inst + r_zp

    du = 99
    dg = 99
    dr = 99

    # Due to the interdependancy of the equations, all bands must be done simultaneously.
    iteration = 0
    while du + dg + dr > 0.0001:
        iteration += 1

        # u calculation
        u_sdss_new = u_inst + u_zp + a_u*(u_sdss - g_sdss)
        # g calculation
        g_sdss_new = g_inst + g_zp + a_g*(g_sdss - r_sdss)
        # r calculation
        r_sdss_new = r_inst + r_zp + a_r*(g_sdss - r_sdss)

        du = abs(u_sdss - u_sdss_new)
        dg = abs(g_sdss - g_sdss_new)
        dr = abs(r_sdss - r_sdss_new)

        u_sdss = u_sdss_new
        g_sdss = g_sdss_new
        r_sdss = r_sdss_new

        print("Iteration {:>03d} | du: {:>06.3f} | dg: {:>06.3f} | dr: {:>06.3f} |".format(iteration, du, dg, dr), end='\n')


    print("\n\nConverged!")
    print("u_sdss: {:.3f}".format(u_sdss))
    print("g_sdss: {:.3f}".format(g_sdss))
    print("i_sdss: {:.3f}".format(i_sdss))

    print("This is a colour term of...")
    print("a_u*(u_sdss - g_sdss) = {:.3f}".format(a_u*(u_sdss - g_sdss)))
    print("a_g*(g_sdss - r_sdss) = {:.3f}".format(a_g*(g_sdss - r_sdss)))
    print("a_i*(g_sdss - r_sdss) = {:.3f}".format(a_i*(g_sdss - r_sdss)))
    print("\n\n\n\n")

    sdss_mags['u'].append(u_sdss)
    sdss_mags['g'].append(g_sdss)
    sdss_mags['r'].append(r_sdss)

sdss_df = pd.DataFrame()
for ap_index in range(len(mags)):
    row = {
        "aperture": str(ap_index+1),
        'u': sdss_mags['u'][ap_index],
        'g': sdss_mags['g'][ap_index],
        'i': sdss_mags['r'][ap_index],
    }
    sdss_df = sdss_df.append(row, ignore_index=True)
print("Comparison star SDSS magnitudes:")
print(sdss_df)

sdss_df.to_csv(oname)
