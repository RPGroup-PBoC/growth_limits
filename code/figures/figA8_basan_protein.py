import numpy as np
import pandas as pd
from scipy import stats
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker
import prot.viz
import prot.size
colors, palette = prot.viz.bokeh_theme()
# dataset_colors = prot.viz.dataset_colors()
prot.viz.plotting_style()

# Exponential fit function
from scipy.optimize import curve_fit
def func(x, a, c, d):
    return a*np.exp(-c*x)+d

basan_df = pd.read_csv('../../data/basan2015_raw_data/basan2015_data.csv')

##############################################
##############################################
# Perform exponential fit of Basan data
##############################################
##############################################

popt_fg, pcov_fg = curve_fit(func, basan_df.growth_rate_hr.values, basan_df.protein_fg.values, p0=(1, 1e-6, 1))

##############################################
##############################################
# Now plot!
##############################################
##############################################

fig, ax = plt.subplots(1, 1, figsize=(2.5, 2))

x = np.linspace(0,2,100)

ax.plot(basan_df.growth_rate_hr.values, basan_df.protein_fg.values, 'o', ms=4, color=colors['light_purple'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'Basan et al 2015')
ax.plot(x, func(x, *popt_fg), '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax.set_ylabel('protein mass per cell [fg]', fontsize=6)
ax.legend(loc = 'upper left', fontsize=6)

ax.xaxis.set_tick_params(labelsize=5)
ax.yaxis.set_tick_params(labelsize=5)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_xlim([0, 2])


plt.tight_layout()
fig.savefig('../../figures/figA8_protein_basan.pdf', bbox_inches='tight')
