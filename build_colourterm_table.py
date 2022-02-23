import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import minimize

instrument = 'ucam'
telescope = 'ntt'
logg = '850'
mydir = os.path.dirname(__file__)

MIST_df = pd.read_csv(os.path.join(mydir, "tables/MIST_mags.csv"))
koester_df = pd.read_csv(os.path.join(mydir, "tables/koester_magnitudes_logg_{}.csv".format(logg)))

diagnostics = ['u-g', 'g-r', 'g-i']
bands = ['u', 'g', 'r', 'i']


# Header
print("{:^12s}&{:^12s}&{:^17s}&{:^17s}\\\\".format("Correction", "Diagnostic", "y-intercept", "Colour Gradient"))
print(r"\hline")

for targetband in bands:
    writeband = "${0}-{0}_s$".format(targetband)
    print(r'\hline')
    for diagnostic in diagnostics:
        def chisq(args):
            '''Uses this model to calculate chi squared:

            g_sdss - g_inst = g_zp + a_g(g-r)_sdss

            optimises for g_zp and a_g
            '''
            chisq = 0.0
            zero_point, colour_term = args

            calc_wd_correction = zero_point + (colour_term * koester_df[diagnostic])
            diff = calc_wd_correction - koester_df['{0}-{0}_s'.format(targetband)]
            chisq += (diff**2).sum()

            color = MIST_df[diagnostic]
            calc_ms_correction = zero_point + (colour_term * color)
            diff = calc_ms_correction - (MIST_df['sdss:{0}'.format(targetband)] - MIST_df['{}:{}:{}_s'.format(instrument, telescope, targetband)])
            chisq += (diff**2).sum()

            return chisq


        # Run the minimiser
        x0 = np.array([24.0, 0.1])
        soln = minimize(
            chisq, x0
        )
        xf = soln['x']
        zero_point, colour_term = xf

        # Table output
        write_diag = "$({})$".format(diagnostic)
        print("{:^12s}&{:^12s}&${:^ 17.3f}$&${:^ 17.3f}$\\\\".format(writeband, write_diag, zero_point, colour_term))

        writeband = ''
print(r'\hline')
print(r'\hline')
