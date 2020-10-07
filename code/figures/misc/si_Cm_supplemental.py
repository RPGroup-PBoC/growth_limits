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
                'size' : d['cell size (μm3)'].mean()}
    data_si_mean = data_si_mean.append(data_list,
                                      ignore_index = True)

# something weird about one data point -. suggests very high number ori; might be due to
# too high Cm and very low growth rate making ht measurement difficult
data_si_mean = data_si_mean[data_si_mean['C+D period (minutes)' ] != 417.6660]

# data_si_mean = data_si_mean[data_si_mean.strain != 'tCRISPRi (MG1655)']

data_si_mean_ = data_si_mean[data_si_mean['type of perturbation'] == 'chloramphenicol']
color_dict = dict(zip(data_si_mean_['growth media'].unique(), palette))
print(data_si_mean_['growth media'].unique())
markers_dict = dict(zip(data_si_mean_['growth media'].unique(),
    ["o", "h", "^", "s", "X", "D", ">", "*", "P", "v", "p"]))
data_si_mean_['concentration'] = data_si_mean_['concentration'].values.astype(np.float)
data_si_mean_ = data_si_mean_.sort_values(by='concentration')



# plot of RNA/protein vs. num ori/num ter
# additional plot of Si chlor data
fig, ax = plt.subplots(2, 2, figsize = (10,8))
ax = ax.ravel()

ax[0].axis('off')


for c, d in data_si_mean_.groupby(['type of perturbation', 'concentration', 'growth media'], sort=False):

    if c[0] != 'chloramphenicol':
        continue

    if c[1] > 7.0:
        k = '#7AA974'
        label = '8-15 $\mu$M'
    if 1.0 < c[1] <= 7.0:
        k = '#D56C55'
        label = '2-7 $\mu$M'
    if c[1] <= 1:
        k = '#EAC264'
        label = '0.15-1 $\mu$M'
    tau_C = d['C period (minutes)' ].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori_ter = 2**(tau_C /tau)
    ax[1].plot(ori_ter,  d['RNA/protein']/2.1, marker = markers_dict[d['growth media'].unique()[0]], color= k,#'k',
                    lw = 0, alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10, label = label)



for c, d in data_si_mean.groupby(['type of perturbation', 'concentration', 'strain'], sort=False):
    if c[0] != 'nutrient conditions':
        continue

    k = '#EAC264'
    tau_C = d['C period (minutes)' ].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori_ter = 2**(tau_C /tau)
    ax[1].plot(ori_ter,  d['RNA/protein']/2.1, 'o', color= k,#'k',
                    alpha=0.5, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10, label = [c[1], c[2]])

ax[1].set_xlabel('# ori/ # ter')
ax[1].set_ylabel('ribosomal fraction\n(RNA/protein ratio x 2.1)')
# ax[1].legend()
##################################
##################################




for c, d in data_si_mean_.groupby(['type of perturbation', 'growth media']):
    if c[0] != 'chloramphenicol':
        continue

    tau_cyc = d['C+D period (minutes)' ].values
    tau_c = d['C period (minutes)' ].values
    # tau = d['doubling time (minutes)'].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori = 2**(tau_cyc /tau)

    # predict total dry mass in fg (1.1 g/ml mass density; 30% dry mass)
    pred_protein = ((1.1 * 0.3*0.55)*(d['size']*1E-12)*1E15)
    Naa = pred_protein * 1E-15 * 6.022E23 / 110
    Naa_ribosomes = Naa * (d['RNA/protein']/2.1)
    N_ribosomes = Naa_ribosomes / 7459.0

    ax[2].plot(ori, N_ribosomes, marker = markers_dict[c[1]], color= color_dict[c[1]],#k,
                    lw = 0, alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10, label = c[1])

for c, d in data_si_mean.groupby(['type of perturbation', 'concentration', 'strain'], sort=False):

    if c[0] != 'nutrient conditions':
        continue

    k = '#EAC264'
    tau_cyc = d['C+D period (minutes)' ].values
    tau_c = d['C period (minutes)' ].values
    # tau = d['doubling time (minutes)'].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori = 2**(tau_cyc /tau)

    # predict total dry mass in fg (1.1 g/ml mass density; 30% dry mass)
    pred_protein = ((1.1 * 0.3*0.55)*(d['size']*1E-12)*1E15)
    Naa = pred_protein * 1E-15 * 6.022E23 / 110
    Naa_ribosomes = Naa * (d['RNA/protein']/2.1)
    N_ribosomes = Naa_ribosomes / 7459.0

    ax[2].plot(ori,  N_ribosomes, 'o', color= k,
                    alpha=0.5, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

ax[2].set_xlim(0,17)
# ax[2].set_ylim(0,90000)

ax[2].set_xlabel('# ori')
ax[2].set_ylabel('estimated ribosome copy number')
ax[2].legend()
###########################################
###########################################

for c, d in data_si_mean_.groupby(['type of perturbation', 'growth media']):
    if c[0] != 'chloramphenicol':
        continue

    tau_cyc = d['C+D period (minutes)' ].values
    tau_c = d['C period (minutes)' ].values
    # tau = d['doubling time (minutes)'].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori = 2**(tau_cyc /tau)

    # predict total dry mass in fg (1.1 g/ml mass density; 30% dry mass)
    pred_protein = ((1.1 * 0.3*0.55)*(d['size']*1E-12)*1E15)
    Naa = pred_protein * 1E-15 * 6.022E23 / 110
    Naa_ribosomes = Naa * (d['RNA/protein']/2.1)
    N_ribosomes = Naa_ribosomes / 7459.0

    ax[3].plot(tau, N_ribosomes/ori, marker = markers_dict[c[1]], color = color_dict[c[1]],
                    lw = 0.75, ls = '--', alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

# for c, d in data_si_mean.groupby(['type of perturbation', 'concentration', 'strain'], sort=False):
#
#     if c[0] != 'nutrient conditions':
#         continue
#
#     k = '#EAC264'
#     tau_cyc = d['C+D period (minutes)' ].values
#     tau_c = d['C period (minutes)' ].values
#     # tau = d['doubling time (minutes)'].values
#     tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
#     ori = 2**(tau_cyc /tau)
#
#     # predict total dry mass in fg (1.1 g/ml mass density; 30% dry mass)
#     pred_protein = ((1.1 * 0.3*0.55)*(d['size']*1E-12)*1E15)
#     Naa = pred_protein * 1E-15 * 6.022E23 / 110
#     Naa_ribosomes = Naa * (d['RNA/protein']/2.1)
#     N_ribosomes = Naa_ribosomes / 7459.0
#
#     ax[3].plot(tau,  N_ribosomes/ori, 'o', color= k,
#                     alpha=0.5, markeredgecolor='k', markeredgewidth=0.25,
#                     ms=4, zorder=10, label = [c[1], c[2]])


ax[3].set_ylabel('estimated number of ribosomes / #ori')
ax[3].set_xlabel('doubling time')


fig.savefig('../../figures/supplemental_Si_Cm.pdf', bbox_inches = 'tight')
