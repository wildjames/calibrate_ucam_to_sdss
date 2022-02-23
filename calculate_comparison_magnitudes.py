import argparse
import json
from os.path import abspath, join, split

import hipercam as hcam
import pandas as pd

from constructReference import get_instrumental_mags

desc = '''
Using the standard - like reduction, calculate the comparison star magnitudes of a logfile.
g_sdss = g_inst,noatmos - g_zp - a_g(g-r)
'''


parser = argparse.ArgumentParser(description=desc)
parser.add_argument("reduction")
# parser.add_argument("RA")
# parser.add_argument("DEC")
parser.add_argument('--oname', default="comparison_star_sdss_mags", required=False)

args = parser.parse_args()

print(args)

fname = args.reduction
target_coords = input("Please enter [RA DEC], space separated, of the target stars: ") # "{} {}".format(args.RA, args.DEC)
oname = args.oname

obsname = 'lasilla'
k_ext = [0.1129, 0.2020, 0.4868]
bands = input("Please enter the observation bands, space separated, in CCD order: ").split(' ') #['r', 'g', 'u']

if 'r' in bands:
    redband = 'r'
elif 'i' in bands:
    redband = 'i'


for i, band in enumerate(bands):
    k_ext[i] = float(input("Enter extinction for {} band: ".format(band)))

obsname = input("Please enter the observation site name: ")

# Gather some data files
# mydir = split(abspath(__file__))[0]
values_fname = abspath("FOUND_VALUES.json")

print("and the prior variables stored at: ")
print(values_fname)

# First, remember what we've found so far
with open(values_fname, 'r') as f:
    variables = json.load(f)

print("I'm using the following variables for colour correction:")
print(variables)


data = hcam.hlog.Hlog.rascii(fname)
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
a_r = variables['a{}'.format(redband)]

# Zero points in each CCD for UCAM on the NTT.
u_zp = variables['u_zp']
g_zp = variables['g_zp']
r_zp = variables['{}_zp'.format(redband)]


for ap_index in range(len(inst_mags)):
    u_inst = inst_mags['u'][ap_index]
    g_inst = inst_mags['g'][ap_index]
    r_inst = inst_mags[redband][ap_index]


    print("Got the following instrumental magnitudes for aperture {}:".format(ap_index+1))
    print("{}: {:.3f}".format(redband, r_inst))
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
    while du + dg + dr > 0.00001:
        iteration += 1

        # u calculation
        u_sdss_new = u_inst + u_zp + a_u*(u_sdss - g_sdss)
        # g calculation
        g_sdss_new = g_inst + g_zp + a_g*(g_sdss - r_sdss)
        # r calculation
        r_sdss_new = r_inst + r_zp + a_r*(g_sdss - r_sdss)

        du = abs(u_sdss - u_sdss_new)
        dg = abs(g_sdss - g_sdss_new)
        dr = abs(r_sdss - r_sdss_new)

        u_sdss = u_sdss_new
        g_sdss = g_sdss_new
        r_sdss = r_sdss_new

        print("Iteration {:>03d} | du: {:>06.3f} | dg: {:>06.3f} | dr: {:>06.3f} |".format(iteration, du, dg, dr), end='\n')


    print("\n\nConverged!")
    print("u_inst: {:.3f}".format(u_inst))
    print("g_inst: {:.3f}".format(g_inst))
    print("{}_inst: {:.3f}".format(redband, r_inst))
    print("")
    print("u_zp:   {:.3f}".format(u_zp))
    print("g_zp:   {:.3f}".format(g_zp))
    print("{}_zp:   {:.3f}".format(redband, r_zp))
    print("")
    print("u_sdss: {:.3f}".format(u_sdss))
    print("g_sdss: {:.3f}".format(g_sdss))
    print("{}_sdss: {:.3f}".format(redband, r_sdss))

    print("This is a colour term of...")
    print("a_u*(u_sdss - g_sdss) = {:.3f}".format(a_u*(u_sdss - g_sdss)))
    print("a_g*(g_sdss - {}_sdss) = {:.3f}".format(redband, a_g*(g_sdss - r_sdss)))
    print("a_r*(g_sdss - {}_sdss) = {:.3f}".format(redband, a_r*(g_sdss - r_sdss)))
    print("\n\n\n\n")

    sdss_mags['u'].append(u_sdss)
    sdss_mags['g'].append(g_sdss)
    sdss_mags[redband].append(r_sdss)

sdss_df = pd.DataFrame()
for ap_index in range(len(inst_mags)):
    row = {
        "aperture": str(ap_index+1),
#       'u': sdss_mags['u'][ap_index],
#       'g': sdss_mags['g'][ap_index],
#       'r': sdss_mags['r'][ap_index],
    }
    for band in bands:
        row[band] = sdss_mags[band][ap_index]
    sdss_df = sdss_df.append(row, ignore_index=True)
print("Comparison star SDSS magnitudes:")
print(sdss_df)

sdss_df.to_csv(oname)
