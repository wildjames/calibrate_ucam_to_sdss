import numpy as np
from calphot.constructReference import get_instrumental_mags
import hipercam as hcam

# Observed instrumental super NTT mags
# Aperture 5
# u_inst = -6.170
# g_inst = -9.131
# r_inst = -9.207
fname = 'Quality_Reductions/2019_09_27-run015.log'

data = hcam.hlog.Hlog.read(fname)
coords = "20 29 17.13 -43 40 19.8"
k_ext = [0.1129, 0.2020, 0.4868]
obsname = 'lasilla'

mags = get_instrumental_mags(data, coords, obsname, k_ext)
u_inst = mags['3'][4]
g_inst = mags['2'][4]
r_inst = mags['1'][4]
print("Got the following instrumental magnitudes:")
print("r: {:.3f}".format(r_inst))
print("g: {:.3f}".format(g_inst))
print("u: {:.3f}".format(u_inst))

# Found parameters. Don't touch unless you know what you're doing!
# Colour terms. a_u is u-g, g and r are g-r.
a_u = -0.037
a_g = -0.024
a_r = -0.032

# Zero points in each CCD for UCAM on the NTT.
u_zp = 24.817
g_zp = 26.218
r_zp = 25.785

# First iteration uses no colour term, so set the sdss mags as equal to the instrumental ones.
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

    print("Iteration {:>03d} | du: {:>06.3f} | dg: {:>06.3f} | dr: {:>06.3f} |".format(iteration, du, dg, dr), end='\r')


print("\n\nConverged!")
print("u_sdss: {:.3f}".format(u_sdss))
print("g_sdss: {:.3f}".format(g_sdss))
print("r_sdss: {:.3f}".format(r_sdss))

print("This is a colour term of...")
print("a_u*(u_sdss - g_sdss) = {:.3f}".format(a_u*(u_sdss - g_sdss)))
print("a_g*(g_sdss - r_sdss) = {:.3f}".format(a_g*(g_sdss - r_sdss)))
print("a_r*(g_sdss - r_sdss) = {:.3f}".format(a_r*(g_sdss - r_sdss)))
