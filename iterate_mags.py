import numpy as np

# Observed instrumental super NTT mags
# SA 114 548 (verification step)
# u_inst = -9.309
# g_inst = -13.936
# r_inst = -14.614

# Aperture 2
u_inst = -4.9831
g_inst = -9.3194
r_inst = -9.9180

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
