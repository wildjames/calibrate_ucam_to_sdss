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

highlight_cells = (
    (0,0),
    (1,1),
    (2,1),
    (3,2),
)



def plot_colterm(ax, diagnostic, targetband):

    # Pre-calculate the MIST colours, since I was too dumb to do it before >:0
    MIST_df['{0}-{0}_s'.format(targetband)] = MIST_df['sdss:{0}'.format(targetband)] - MIST_df['{}:{}:{}_s'.format(instrument, telescope, targetband)]

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
        diff = calc_ms_correction - (MIST_df['{0}-{0}_s'.format(targetband)])
        chisq += (diff**2).sum()

        return chisq

    # Run the minimiser
    x0 = np.array([24.0, 0.1])
    soln = minimize(
        chisq, x0
    )
    print("Fitting a straight line resulted in the following:")
    print(soln)


    MIST_df.plot.scatter(
        diagnostic, '{0}-{0}_s'.format(targetband),
        ax=ax,
        color='black',
        s=3
        #label="MIST MS models, 8.5Gyrs",
    )
    koester_df.plot.scatter(
        diagnostic, '{0}-{0}_s'.format(targetband),
        ax=ax,
        color='red',
        s=3
        #label='WD Koester Models, logg {}'.format(logg),
    )

    # Best fit line
    zp, colterm = soln['x']

    xr = np.linspace(MIST_df[diagnostic].min(), MIST_df[diagnostic].max(), 100)
    yr = zp + colterm * xr

    ax.plot(xr, yr, color='blue') #label='Fitted line, colour term: {:.3f}'.format(colterm))

    return ax

labels = {
    "u-g": "u'-g'",
    "g-r": "g'-r'",
    "g-i": "g'-i'",
}

ny, nx = len(bands), len(diagnostics)
fig, ax = plt.subplots(ny,nx, figsize=(8,6))
for i, diag in enumerate(diagnostics):
    for j, targ in enumerate(bands):
        ax[j, i] = plot_colterm(ax[j,i], diag, targ)
        ax[j, i].yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        ax[j, i].xaxis.set_major_formatter(FormatStrFormatter("%.0f"))

        label = labels[diag]
        ax[j, i].set_xlabel("$({})$".format(label))
        ax[j, i].set_ylabel("$({0}'-{0}".format(targ)+r"_{sup})$")

        if j!=3:
            ax[j, i].set_xlabel('')
            ax[j, i].set_xticks([], [])
        if i!=0:
            ax[j, i].set_ylabel('')
            ax[j, i].set_yticks([], [])

for i,j in highlight_cells:
    ax[i,j].set_facecolor('lightgrey')

fig.subplots_adjust(right=0.94, top=0.94, hspace=0, wspace=0)
plt.savefig("figs/colour_term_tracks.pdf")
plt.show()
plt.close()
