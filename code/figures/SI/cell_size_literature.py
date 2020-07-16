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

# Exponential fit functions
from scipy.optimize import curve_fit
def func(x, a, c, d):
    return a*np.exp(-c*x)+d


##############################################
##############################################
# gather cell size data
##############################################
##############################################

size_df = pd.read_csv('../../../data/cell_size_literature/size_data.csv')

##############################################
##############################################
# Perform exponential fit of Basan data for reference
#  (we have fit curves for Si 2017 and T-A 2015 from
#  their papers)
##############################################
##############################################

popt_vol, pcov_vol = curve_fit(func, size_df[size_df.dataset_name == 'Basan et al 2015'].growth_rate_hr.values,
                size_df[size_df.dataset_name == 'Basan et al 2015'].volume_um3.values, p0=(1, 1e-6, 1))

##############################################
##############################################
# Now plot!
##############################################
##############################################
#
# fig, ax = plt.subplots(1, 2, figsize=(5, 2.5))

fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))

colordict = {'Si et al 2017' : colors['pale_red'],
            'Taheri-Araghi et al 2015' : colors['pale_green'],
            'Basan et al 2015' : colors['light_purple'],
            'Schmidt et al 2016' : colors['light_blue']}

x = np.linspace(0,2.5,150)
for ax_ in ax:
    for c, d in size_df.groupby(['dataset_name']):
        ax_.plot(d['growth_rate_hr'].values, d['volume_um3'].values, 'o', ms=4, color=colordict[c],
                markeredgewidth=0.5, markeredgecolor='k', label = c)
        ax_.set_ylabel(r'average cell size [$\mu$m$^3$]', fontsize=6)

        ax_.plot(x, 0.28 * np.exp(1.33  * x), '-', alpha = 0.6,
                    lw = 0.5, color = 'k', zorder = 0)
        ax_.plot(x, 0.27 * 2**(1.1  * x / np.log(2)), '-', alpha = 0.6,
                    lw = 0.5, color = 'k', zorder = 0)
        ax_.plot(x[:-30], func(x[:-30], *popt_vol), '-', alpha = 0.6,
                    lw = 0.5, color = 'k', zorder = 0)

for ax_ in ax:
    ax_.xaxis.set_tick_params(labelsize=5)
    ax_.yaxis.set_tick_params(labelsize=5)
    ax_.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    ax_.set_xlim([0, 2.5])

ax[1].set_yscale('log')
ax[1].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

ax[1].set_yticks([0.3, 0.5, 1, 2, 3, 5])

handles, labels = ax[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax[0].legend(by_label.values(), by_label.keys(), loc = 'upper left', fontsize = 6)


plt.tight_layout()
fig.savefig('../../../figures/size_data_lit.pdf', bbox_inches='tight')
