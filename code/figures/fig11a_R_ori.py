#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import prot.viz
import prot.size as size
colors, palette = prot.viz.bokeh_theme()

dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
prot.viz.plotting_style()

from scipy.optimize import curve_fit
# def func(x, a, c, d):
#     return a*np.exp(-c*x)+d

######################
# grab data from Si 2017 to extimate <# ori>
######################
# Load the full Si 2017 SI
data_si = pd.read_csv('../../data/si_2017_raw/si2017_full.csv')

# consider the mean values for each strain/condition they considered
data_si_mean = pd.DataFrame()
for c, d in data_si.groupby(['strain type (background)', 'growth media', 'inducer or drug added', 'concentration ', 'type of perturbation']):
    data_list = {'strain' : c[0],
                 'growth media' : c[1],
                 'inducer or drug added' : c[2],
                 'concentration': c[3],
                 'number of origins' :d['number of origins'].mean(),
                 'RNA/protein' : d['RNA/protein'].mean(),
                 'DNA content per cell (genome equivalents)' : d['DNA content per cell (genome equivalents)'].mean(),
                 'type of perturbation' :c[4],
                 'C+D period (minutes)' : d['C+D period (minutes)'].mean(),
                 'C period (minutes)' : d['C period (minutes)' ].mean(),
                 'growth_rate_hr' : d['growth rate (1/hours)'].mean(),
                 'doubling time (minutes)' : d['doubling time (minutes)'].mean() }
    data_si_mean = data_si_mean.append(data_list,
                                      ignore_index = True)

# perform fit of growth rate vs. t_cyc so that we can estimate number of
# origins (2**(t_cyc/tau)) where tau is the doubling time.
data_si_mean = data_si_mean.sort_values(by='growth_rate_hr')
data_si_mean['tau'] = 60*(np.log(2)/data_si_mean['growth_rate_hr'])
data_si_mean = data_si_mean[data_si_mean.strain != 'tCRISPRi (MG1655)']
data_si_mean = data_si_mean[data_si_mean['type of perturbation'] == 'nutrient conditions']


def func_lin(l, a, b):
    x = 60*(np.log(2)/l)
    return a*x + b

#### piecewise fit for t_C- piecewise fit works well with transition at 40 min
t_C_const = data_si_mean[data_si_mean.tau <=40]['C period (minutes)'].mean()
t_C_lin = data_si_mean[data_si_mean.tau > 40]['C period (minutes)'].values
l_lin =  data_si_mean[data_si_mean.tau > 40]['growth_rate_hr'].values
popt_tC_lin, pcov_tC_lin = curve_fit(func_lin, l_lin, t_C_lin, p0=(1, 1))

#### piecewise fit for t_cyc - piecewise fit works well with transition at 43 min
t_cyc_const = data_si_mean[data_si_mean.tau <=43]['C+D period (minutes)'].mean()
t_cyc_lin = data_si_mean[data_si_mean.tau > 43]['C+D period (minutes)'].values
l_lin =  data_si_mean[data_si_mean.tau > 43]['growth_rate_hr'].values
popt_tcyc_lin, pcov_tcyc_lin = curve_fit(func_lin, l_lin, t_cyc_lin, p0=(1, 1))


# %%
######################
# plot 2  #
######################
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# plot of RNA/protein vs. num ori
# Load the complex subunit counts.
subunits = pd.read_csv('../../data/compiled_annotated_complexes.csv')

# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# Compute the minimum number of complexes.
complex_count = subunits.groupby(['dataset', 'dataset_name', 'condition',
                    'growth_rate_hr', 'complex_annotation',
                    'complex'])['n_units'].mean().reset_index()

complex_ribo = complex_count[complex_count.complex_annotation == 'ribosome']

complex_ribo['tau']  = 60*(np.log(2)/(complex_ribo['growth_rate_hr']))
t_cyc_arr = []
for i, val in enumerate(complex_ribo['tau'].values):
    if val <= 40:
        t_cyc_ = t_cyc_const
    else:
        t_cyc_ = func_lin(complex_ribo['growth_rate_hr'].values[i], *popt_tcyc_lin)
    t_cyc_arr = np.append(t_cyc_arr, t_cyc_)

complex_ribo['t_cyc'] = t_cyc_arr #func_lin(complex_ribo['growth_rate_hr'], *popt_tcyc_lin)
complex_ribo['# ori'] =  2**(complex_ribo['t_cyc'] / complex_ribo['tau'] )

# for g, d in complex_ribo.groupby(['dataset', 'condition', 'growth_rate_hr']):
for g, d in complex_ribo.groupby(['dataset', 'dataset_name']):
    ax.plot(d['# ori'], d['n_units'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[1], ms=4, zorder=10)

ax.set_xlabel('estimated # ori', fontsize=6)
ax.set_ylabel('ribosomes per cell', fontsize=6)
ax.xaxis.set_tick_params(labelsize=5)
ax.yaxis.set_tick_params(labelsize=5)
ax.legend(fontsize=6, loc = 'upper left')
plt.savefig('../../figures/fig11a_ribosome_vs_ori.pdf', bbox_inches='tight')
