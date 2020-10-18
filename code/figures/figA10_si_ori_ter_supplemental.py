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
def func_lin(l, a, b):
    x = 60*(np.log(2)/l)
    return a*x + b

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


#### piecewise fit for t_C- piecewise fit works well with transition at 40 min
t_C_const = data_si_mean[data_si_mean.tau <=40]['C period (minutes)'].mean()
t_C_lin = data_si_mean[data_si_mean.tau > 40]['C period (minutes)'].values
l_lin =  data_si_mean[data_si_mean.tau > 40]['growth_rate_hr'].values
popt_tC_lin, pcov_tC_lin = curve_fit(func_lin, l_lin, t_C_lin, p0=(1, 1))
print('t_C: constant regime, ' , t_C_const, 'linear fit:', popt_tC_lin)

#### piecewise fit for t_cyc - piecewise fit works well with transition at 43 min
t_cyc_const = data_si_mean[data_si_mean.tau <=43]['C+D period (minutes)'].mean()
t_cyc_lin = data_si_mean[data_si_mean.tau > 43]['C+D period (minutes)'].values
l_lin =  data_si_mean[data_si_mean.tau > 43]['growth_rate_hr'].values
popt_tcyc_lin, pcov_tcyc_lin = curve_fit(func_lin, l_lin, t_cyc_lin, p0=(1, 1))
print('t_cyc: constant regime, ' , t_cyc_const, 'linear fit:', popt_tcyc_lin)

# %%
######################
# plot configuration #
######################
fig = plt.figure(figsize = (6,4))
widths = [2, 2]
heights = [2, 2]
spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
                          height_ratios=heights)

ax1 = fig.add_subplot(spec[0])
ax2 = fig.add_subplot(spec[1])
ax3 = fig.add_subplot(spec[2])
ax4 = fig.add_subplot(spec[3])


for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'nutrient conditions':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    ax1.plot(d['tau'],  d['C period (minutes)'], 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)
    ax2.plot(d['tau'],  d['C+D period (minutes)'], 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

    ax3.plot(d['growth_rate_hr'],  d['C period (minutes)'], 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)
    ax4.plot(d['growth_rate_hr'],  d['C+D period (minutes)'], 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

# plot 1 , 3; t_C vs. tau and vs. lambda
tau = np.linspace(0.1, 40, 100)
x = (np.log(2)/tau)*60
tC_x = t_C_const.mean()*np.ones(100)
ax1.plot(tau,  tC_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')
ax3.plot(x,  tC_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')

tau = np.linspace(40, (np.log(2)/0.1)*60, 100)
x = (np.log(2)/tau)*60
tC_x = func_lin(x, *popt_tC_lin)
ax1.plot(tau,  tC_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')
ax3.plot(x,  tC_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')

# plot 2 ; t_cyc vs. tau and lambda
tau = np.linspace(0.1, 40, 100)
x = (np.log(2)/tau)*60
tcyc_x = t_cyc_const.mean()*np.ones(100)
ax2.plot(tau,  tcyc_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')
ax4.plot(x,  tcyc_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')


tau = np.linspace(40, (np.log(2)/0.1)*60, 100)
x = (np.log(2)/tau)*60
tcyc_x = func_lin(x, *popt_tcyc_lin)
ax2.plot(tau,  tcyc_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')
ax4.plot(x,  tcyc_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')

for a in [ax1,ax2]:
    a.set_xlabel('doubling time [min]', fontsize=6)
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
    a.set_xlim([0, 130])
    a.set_ylim([0, 150])

for a in [ax3,ax4]:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
    a.set_xlim([0, 2])
    a.set_ylim([0, 150])

ax1.set_ylabel('t$_{C}$ [min]', fontsize=6)
ax3.set_ylabel('t$_{C}$ [min]', fontsize=6)
ax2.set_ylabel('t$_{cyc}$ [min]', fontsize=6)
ax4.set_ylabel('t$_{cyc}$ [min]', fontsize=6)

fig.savefig('../../figures/figA10_supplemental_ori_ter.pdf', bbox_inches = 'tight')
