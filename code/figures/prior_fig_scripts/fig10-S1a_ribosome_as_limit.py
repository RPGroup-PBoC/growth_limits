#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

from scipy.optimize import curve_fit
def func(x, a, c, d):
    return a*np.exp(-c*x)+d

# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

L_R = 7459.0 # length of all subunits in ribosomes, in amino acids

ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']

# %%
######################
# plot configuration #
######################
fig = plt.figure(figsize = (6,5))
widths = [5, 2, 2]
heights = [2, 4, 6]
spec = fig.add_gridspec(ncols=3, nrows=3, width_ratios=widths,
                          height_ratios=heights)

ax1 = fig.add_subplot(spec[2, 0])

# %%
######################
# ribosomal fraction information#
######################
# determine  ribosome fraction in proteomic data sets
df_ribo_frac = pd.DataFrame()
for c, d in data.groupby(['dataset', 'condition', 'growth_rate_hr', 'dataset_name']):
    mass_ribo = d[d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    frac_ribo = (mass_ribo )/ d.fg_per_cell.sum()

    data_list = {'frac_ribo' : frac_ribo,
                'dataset' : c[0],
                'condition' : c[1],
                'growth_rate_hr' : c[2],
                'dataset_name' : c[3]}

    df_ribo_frac = df_ribo_frac.append(data_list,
                                        ignore_index = True)

# add in data from Scott et al.
[lambda_scott, R_P_scott] = np.array(\
        [[0.4, 0.177],
        [0.57, 0.230],
        [0.71, 0.221],
        [1.00, 0.287],
        [1.31, 0.414],
        [1.58, 0.466]]).T

# Add in data from Forchhammer & Lindahl
[lambda_for, R_P_for] = np.array(\
        [[0.38, 0.189],
        [0.60, 0.224],
        [1.04, 0.295],
        [1.46, 0.421],
        [1.73, 0.469]]).T

# Bremmer and Dennis
[lambda_brem, R_P_brem] = np.array(\
        [[0.42, 0.200],
        [0.69, 0.225],
        [1.04, 0.331],
        [1.39, 0.391],
        [1.73, 0.471]]).T

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
                 'growth_rate_hr' : d['growth rate (1/hours)'].mean()}
    data_si_mean = data_si_mean.append(data_list,
                                      ignore_index = True)

# load in the Dai et al. 2016 data
dai_nut_df = pd.read_csv('../../data/dai2016_raw_data/dai2016_summary.csv')
dai_nut_df = dai_nut_df[dai_nut_df['Cm (μM)'] == 0]
lambda_dai = dai_nut_df.growth_rate_hr.values
R_P_dai = dai_nut_df.RNA_P_ratio.values

# fit measurements of active fraction from Dai et al. data
dai_nut_df = dai_nut_df.sort_values(by='growth_rate_hr', ascending = True)

popt_dai, pcov_dai = curve_fit(func, dai_nut_df.growth_rate_hr.values, dai_nut_df.f_a.values, p0=(1, 1e-6, 1))

# %%
######################
# Plotting  #
######################

# Add Dai et al, Scott, and historical
ax1.plot((R_P_dai/2.1)*func(lambda_dai, *popt_dai),lambda_dai, 'o', color=  colors['light_yellow'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=0, label = 'Dai et al. 2016')

ax1.plot((R_P_scott/2.1)*func(lambda_scott, *popt_dai),lambda_scott, 'o', color= colors['light_purple'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=0, label = 'Scott et al. 2010')

ax1.plot((R_P_for/2.1)*func(lambda_for, *popt_dai),lambda_for, 'o', color= '#1F4B99',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=0, label = 'Forchhammer & Lindahl 1971')

ax1.plot((R_P_brem/2.1)*func(lambda_brem, *popt_dai),lambda_brem, 'o', color= '#B7741A',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=0, label = 'Bremmer & Dennis 2006')


for g, d in df_ribo_frac.groupby(['dataset', 'condition', 'growth_rate_hr', 'dataset_name']):
    ax1.plot(d['frac_ribo']*func( d.growth_rate_hr.unique(), *popt_dai), g[2], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[3], ms=4, zorder=10)

for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'nutrient conditions':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    ax1.plot((d['RNA/protein']/2.1)*func( d.growth_rate_hr.unique(), *popt_dai),d.growth_rate_hr.unique(), 'o', color= k,
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=0, label = 'Si et al. 2017')


ax1.xaxis.set_tick_params(labelsize=5)
ax1.yaxis.set_tick_params(labelsize=5)

ax1.set_xlim([0, 0.26])
ax1.set_ylim([0, 2.2])
ax1.set_xlabel('estimated active ribosomal\nfraction ($\Phi_R f_a$)', fontsize=6)
ax1.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)

# Plot the prediction.
frac = np.linspace(0,0.3,100)
gr = (17.1 * frac/ L_R) * 3600
ax1.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)

handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc = 'upper left', fontsize = 6)

plt.tight_layout()
plt.savefig('../../figures/fig10-S1a_ribosome_as_limit.pdf')
