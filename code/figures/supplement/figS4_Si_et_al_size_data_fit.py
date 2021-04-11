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

def func2(x, a, c):
    return a*np.exp(-c*x)

##############################################
##############################################
# gather cell size information from Si et al.
# 2017, 2019
##############################################
##############################################

# %% Si et al. 2019 data
Si_filelist = glob.glob('../../data/si_2019_raw/*.csv')

Si_2019 = pd.DataFrame()
for f in Si_filelist:
    Si_temp = pd.read_csv(f)
    Si_temp['condition'] = f.split('supplemental_')[-1].split('.')[0]
    Si_temp['strain type (background)'] = f.split('supplemental_')[-1].split('.')[0].split('_')[0]
    Si_2019 = Si_2019.append(Si_temp)

Si_2019['length'] = Si_2019['newborn size (micron)'] * (2**(Si_2019['elongation rate (1/hour)']*Si_2019['generation time (minute)']/60.0) - 1)/(Si_2019['elongation rate (1/hour)']*(Si_2019['generation time (minute)']/60.0)*np.log(2))
Si_2019['growth_rate_hr'] = np.log(2) * (Si_2019['elongation rate (1/hour)'])

# Size following equation from Basan et al. 2015
Si_2019['size'] = np.pi * (Si_2019['cell width (micron)']/2)**2 * (Si_2019['length'] - 2*(Si_2019['cell width (micron)']/2)/3)

Si_2019_avg = pd.DataFrame()

# calculate average size values for each growth condition/ cell strain
for c, d in Si_2019.groupby('condition'):
    data_list = {'cell size (μm3)' : d['size'].mean(),
                 'cell width (μm)' : d['cell width (micron)'].mean(),
                 'cell length (μm)' : d['length'].mean(),
                 'growth rate (1/hours)' : d.growth_rate_hr.mean(),
                'dataset': 'Si et al. 2019',
                'strain type (background)' : d['strain type (background)'].unique()[0],
                'condition' : c}
    Si_2019_avg = Si_2019_avg.append(data_list,
                                    ignore_index = True)


# %% Si et al. 2017 data
Si_filelist = glob.glob('../../data/si_2017_raw/*.csv')
Si_2017 = pd.DataFrame()
for f in Si_filelist:
    Si_temp = pd.read_csv(f)
    Si_2017 = Si_2017.append(Si_temp)

Si_2017 = Si_2017[Si_2017['type of perturbation'] == 'nutrient conditions']
Si_2017 = Si_2017[Si_2017['strain type (background)'] != 'tCRISPRi (MG1655)']

Si_avg_vol = Si_2019_avg
for c, d in Si_2017.groupby(['experiment name', 'growth media']):
    data_list = {'growth rate (1/hours)' : d['growth rate (1/hours)'].unique()[0],
                 'cell size (μm3)' : d['cell size (μm3)'].unique()[0],
                 'cell length (μm)' : d['cel length (μm)'].unique()[0],
                 'cell width (μm)' : d['cel width (μm)'].unique()[0],
                 'strain type (background)' : d['strain type (background)'].unique()[0],
                 'dataset' : 'Si et al. 2017',
                 'condition' : c[1],
                'experiment' : d['experiment name'].unique()[0]}
    Si_avg_vol = Si_avg_vol.append(data_list,
                         ignore_index = True)



##############################################
##############################################
# Fit strain MG1655 to exponential function
##############################################
##############################################

Si_avg_vol_MG1655 = Si_avg_vol[Si_avg_vol['strain type (background)'] == 'MG1655']

popt_vol_MG1655, pcov_vol_MG1655 = \
        curve_fit(func2, Si_avg_vol_MG1655['growth rate (1/hours)'].values,
                  Si_avg_vol_MG1655['cell size (μm3)'].values, p0=(1, 1e-6))

popt_length_MG1655, pcov_length_MG1655 = \
        curve_fit(func, Si_avg_vol_MG1655[Si_avg_vol_MG1655.experiment !='MG_20160101']['growth rate (1/hours)'].values,
                  Si_avg_vol_MG1655[Si_avg_vol_MG1655.experiment !='MG_20160101']['cell length (μm)'].values, p0=(1, 1e-6,1))

popt_width_MG1655, pcov_width_MG1655 = \
        curve_fit(func2, Si_avg_vol_MG1655[Si_avg_vol_MG1655.experiment !='MG_20160101']['growth rate (1/hours)'].values,
                  Si_avg_vol_MG1655[Si_avg_vol_MG1655.experiment !='MG_20160101']['cell width (μm)'].values, p0=(1, 1e-6))



x = np.linspace(0,2.2,100)
si_avg_fit_MG1655 = pd.DataFrame({'growth rate (1/hours)' : x,
                        'cell size (μm3)' : func2(x, *popt_vol_MG1655),
                        'cell length (μm)' : func(x, *popt_length_MG1655),
                        'cell width (μm)' : func2(x, *popt_width_MG1655),
                        'dataset' : ['fit to Si et al. 2017, 2019; MG1655' for i in np.arange(len(x))]},
                     columns = ['growth rate (1/hours)', 'cell size (μm3)', 'cell length (μm)', 'cell width (μm)',
                               'dataset'])

##############################################
##############################################
# Now plot!
##############################################
##############################################

fig, ax = plt.subplots(1, 3, figsize=(6, 2.5))

marker = {'Si et al. 2019' : 's',
            'Si et al. 2017' : 'o'}

color = {'MG1655' : '#AB85AC',
            'NCM3722' : '#A9BFE3'}


for c, d in Si_avg_vol.groupby(['dataset', 'condition', 'strain type (background)', 'growth rate (1/hours)']):
    # length
    ax[0].plot(c[3], d['cell length (μm)'].values[0], marker[c[0]], ms=4, color=color[c[2]],
            markeredgewidth=0.5, markeredgecolor='k')
    ax[0].set_ylabel(r'average cell length [$\mu$m]', fontsize=6)

    # width
    ax[1].plot(c[3], d['cell width (μm)'].values[0], marker[c[0]], ms=4, color=color[c[2]],
            markeredgewidth=0.5, markeredgecolor='k', label = c[0])
    ax[1].set_ylabel(r'average cell width [$\mu$m]', fontsize=6)

    # size/ volume
    ax[2].plot(c[3], d['cell size (μm3)'].values[0], marker[c[0]], ms=4, color=color[c[2]],
            markeredgewidth=0.5, markeredgecolor='k')
    ax[2].set_ylabel(r'average cell size [$\mu$m$^3$]', fontsize=6)

# fit lines
ax[0].plot(si_avg_fit_MG1655['growth rate (1/hours)'].values, si_avg_fit_MG1655['cell length (μm)'].values,
        color = 'k', lw = 0.5, zorder = 0)
ax[1].plot(si_avg_fit_MG1655['growth rate (1/hours)'].values, si_avg_fit_MG1655['cell width (μm)'].values,
        color = 'k', lw = 0.5, zorder = 0, label = 'fit to strain MG1655')
ax[2].plot(si_avg_fit_MG1655['growth rate (1/hours)'].values, si_avg_fit_MG1655['cell size (μm3)'].values,
        color = 'k', lw = 0.5, zorder = 0)

for ax_ in ax:
    ax_.xaxis.set_tick_params(labelsize=5)
    ax_.yaxis.set_tick_params(labelsize=5)
    ax_.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    ax_.set_yscale('log')
    ax_.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax_.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

    xmin = si_avg_fit_MG1655['growth rate (1/hours)'].values[0]
    xmax = si_avg_fit_MG1655['growth rate (1/hours)'].values[-1]
    ax_.set_xlim([xmin, xmax])

ax[0].set_yticks([2, 3, 4, 5, 7])
ax[1].set_yticks([0.3, 0.5, 1, 1.5])
ax[2].set_yticks([0.3, 0.5, 1, 2, 3, 5])

# add lengend to width plot

# handles, labels = ax[1].get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# legend1 = ax[1].legend(by_label.values(), by_label.keys(), loc = 'lower center', fontsize = 6)


# In total 3x3 lines have been plotted
lines = ax[1].get_lines()

from matplotlib.lines import Line2D

legend_elements1 = [Line2D([0], [0], marker='s', color='gray', label='Si et al. 2019',
                     ms=4, lw = 0, markeredgewidth=0.5, markeredgecolor='k',
                     alpha = 0.4),
                   Line2D([0], [0], marker='o', color='gray', label='Si et al. 2017',
                     ms=4, lw = 0,markeredgewidth=0.5, markeredgecolor='k'),
                   Line2D([0], [0], marker=None, color='k', label='Si et al. 2017',
                     ms=4, lw = 0.5,markeredgewidth=0.5, markeredgecolor='k')]

legend_elements2 = [Line2D([0], [0], marker='o', color='#AB85AC', label='MG1655',
                     ms=4, lw = 0, markeredgewidth=0.5, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='#A9BFE3', label='NCM3722',
                     ms=4, lw = 0,markeredgewidth=0.5, markeredgecolor='k')]

legend1 = ax[1].legend(handles = legend_elements1, loc = 'lower center',
                        fontsize = 6)
legend2 = ax[1].legend(handles = legend_elements2, loc = 'upper left',
                        fontsize = 6)
ax[1].add_artist(legend1)
ax[1].add_artist(legend2)



plt.tight_layout()
fig.savefig('../../figures/figS4_Si_size_data_fit.pdf', bbox_inches='tight')
