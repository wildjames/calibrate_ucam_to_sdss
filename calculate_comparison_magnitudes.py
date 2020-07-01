import hipercam as hcam
import json
from calphot.constructReference import get_instrumental_mags
import pandas as pd


'''
Using the standard - like reduction, calculate the comparison star magnitudes of a logfile.
g_sdss = g_inst,noatmos - g_zp - a_g(g-r) 
'''

fname = "Quality_Reductions/2019_09_27-run015_std.log"
target_coords = "20 29 17.13 -43 40 19.8"
obsname = 'lasilla'
k_ext = [0.1129, 0.2020, 0.4868]
bands = ['r', 'g', 'u']

# Previously calculated values
variables_fname = "FOUND_VALUES.json"
with open(variables_fname, 'r') as f:
    variables = json.load(f)


data = hcam.hlog.Hlog.read(fname)
# Gets a dict of lists, of all the aperture mags. The first one will be the target, so can be discarded.
inst_mags = get_instrumental_mags(data, target_coords, obsname, k_ext)
inst_mags = {bands[int(key)-1]: mags[1:] for key, mags in inst_mags.items()}

print("The 0th aperture is assumed to be the target, so it'll be discarded now.")
print("I got the following comaprison instrumental magnitudes:")
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
a_r = variables['ar']

# Zero points in each CCD for UCAM on the NTT.
u_zp = variables['u_zp']
g_zp = variables['g_zp']
r_zp = variables['r_zp']

for ap_index in range(len(mags)):
    u_inst = inst_mags['u'][ap_index]
    g_inst = inst_mags['g'][ap_index]
    r_inst = inst_mags['r'][ap_index]
    print("Got the following instrumental magnitudes for aperture {}:".format(ap_index+1))
    print("r: {:.3f}".format(r_inst))
    print("g: {:.3f}".format(g_inst))
    print("u: {:.3f}".format(u_inst))

    # First iteration uses no colour term, so set the sdss mags as equal to the instrumental ones.
    u_sdss = u_inst - u_zp
    g_sdss = g_inst - g_zp
    r_sdss = r_inst - r_zp

    du = 99
    dg = 99
    dr = 99

    # Due to the interdependancy of the equations, all bands must be done simultaneously.
    iteration = 0
    while du + dg + dr > 0.0001 or iteration < 10:
        iteration += 1

        # u calculation
        u_sdss_new = u_inst - u_zp - a_u*(u_sdss - g_sdss)
        # g calculation
        g_sdss_new = g_inst - g_zp - a_g*(g_sdss - r_sdss)
        # r calculation
        r_sdss_new = r_inst - r_zp - a_r*(g_sdss - r_sdss)

        du = abs(u_sdss - u_sdss_new)
        dg = abs(g_sdss - g_sdss_new)
        dr = abs(r_sdss - r_sdss_new)

        u_sdss = u_sdss_new
        g_sdss = g_sdss_new
        r_sdss = r_sdss_new

        print("Iteration {:>03d} | du: {:>06.3f} | dg: {:>06.3f} | dr: {:>06.3f} |".format(iteration, du, dg, dr), end='\r')


    print("\n\nConverged!")
    print("u_sdss: {:.3f}".format(u_sdss))
    print("g_sdss: {:.3f}".format(g_sdss))
    print("r_sdss: {:.3f}".format(r_sdss))

    print("This is a colour term of...")
    print("a_u*(u_sdss - g_sdss) = {:.3f}".format(a_u*(u_sdss - g_sdss)))
    print("a_g*(g_sdss - r_sdss) = {:.3f}".format(a_g*(g_sdss - r_sdss)))
    print("a_r*(g_sdss - r_sdss) = {:.3f}".format(a_r*(g_sdss - r_sdss)))
    print("\n\n\n\n")

    sdss_mags['u'].append(u_sdss)
    sdss_mags['g'].append(g_sdss)
    sdss_mags['r'].append(r_sdss)

print(sdss_mags)
sdss_df = pd.DataFrame()
for ap_index in range(len(mags)):
    row = {
        "aperture": str(ap_index+1),
        'u': sdss_mags['u'][ap_index],
        'g': sdss_mags['g'][ap_index],
        'r': sdss_mags['r'][ap_index],
    }
    sdss_df = sdss_df.append(row, ignore_index=True)
print(sdss_df)


