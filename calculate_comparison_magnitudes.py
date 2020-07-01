import hipercam as hcam
import json
from calphot.constructReference import get_instrumental_mags


'''
Using the standard - like reduction, calculate the comparison star magnitudes of a logfile.
g_sdss = g_inst,noatmos - g_zp - a_g(g-r) 
'''

fname = "Quality_Reductions/2019_09_27-run015_std.log"
target_coords = "20 29 17.13 -43 40 19.8"
obsname = 'lasilla'
k_ext = [0.1129, 0.2020, 0.4868]
bands = ['r', 'g', 'u']


data = hcam.hlog.Hlog.read(fname)
# Gets a dict of lists, of all the aperture mags. The first one will be the target, so can be discarded.
inst_mags = get_instrumental_mags(data, target_coords, obsname, k_ext)
inst_mags = {bands[int(key)-1]: mags[1:] for key, mags in inst_mags.items()}

print("I got the following instrumental magnitudes:")
for key, mags in inst_mags.items():
    print("Band: {}".format(key))
    for mag in mags:
        print("  {:.3f}".format(mag))
    print()
