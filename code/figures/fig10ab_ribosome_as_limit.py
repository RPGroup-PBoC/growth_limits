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
fig = plt.figure(figsize = (6,5))#constrained_layout=True)
# widths = [6, 2.5, 2.5, 5]
widths = [5, 2, 2]
heights = [2, 4, 6]
spec = fig.add_gridspec(ncols=3, nrows=3, width_ratios=widths,
                          height_ratios=heights)

# subplot [0,0] blank
# plot of growth rate vs. ribosomal fraction
ax2 = fig.add_subplot(spec[1, 1])

# plot of growth rate vs. active ribosomal fraction
ax3 = fig.add_subplot(spec[2, 0])
# plot of number of ribosomes vs. growth rate - rRNA persepective
ax4 = fig.add_subplot(spec[2, 1:])

# %
############################
# Plot 2 - growth rate vs. ribosomal fraction
############################
#
# # Add in data from Dai et al.

# load in the Dai et al. 2016 data
dai_nut_df = pd.read_csv('../../data/dai2016_raw_data/dai2016_summary.csv')
dai_nut_df = dai_nut_df[dai_nut_df['Cm (μM)'] == 0]
lambda_dai = dai_nut_df.growth_rate_hr.values
R_P_dai = dai_nut_df.RNA_P_ratio.values

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

df_ribo_frac = pd.DataFrame()
for c, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    mass_ribo = d[d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    frac_ribo = (mass_ribo )/ d.fg_per_cell.sum()

    data_list = {'frac_ribo' : frac_ribo,
                'dataset' : c[0],
                'condition' : c[1],
                'growth_rate_hr' : c[2]}

    df_ribo_frac = df_ribo_frac.append(data_list,
                                        ignore_index = True)

#
# Plot the prediction.
frac = np.linspace(0,1.0,100)
gr = (17.1 * frac/ L_R) * 3600
ax2.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)
ax2.hlines((17.1 / L_R) * 3600, 0, np.max(gr), color='k', linestyle='--', lw=0.75, label='__nolegend__')


ax2.xaxis.set_tick_params(labelsize=5)
ax2.yaxis.set_tick_params(labelsize=5)
ax2.set_xlim([0, 1.0])
ax2.set_ylim([0, 9])
ax2.set_xlabel('ribosomal fraction ($\Phi_R$)', fontsize=6)
ax2.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)



# %
############################
# Plot 3 - active ribosomal fraction
############################

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

# fit measurements of active fraction from Dai et al. data
dai_nut_df = dai_nut_df.sort_values(by='growth_rate_hr', ascending = True)

popt_dai, pcov_dai = curve_fit(func, dai_nut_df.growth_rate_hr.values, dai_nut_df.f_a.values, p0=(1, 1e-6, 1))


# Add Dai et al, Scott, and historical
ax3.plot((R_P_dai/2.1)*func(lambda_dai, *popt_dai),lambda_dai, 'o', color=  'k', #color=colors['light_yellow'],
                alpha=0.2, markeredgecolor='k', markeredgewidth=0,
                ms=4, zorder=0, label = 'non-proteomic data')

ax3.plot((R_P_scott/2.1)*func(lambda_scott, *popt_dai),lambda_scott, 'o', color=  'k', #color=colors['light_purple'],
                alpha=0.2, markeredgecolor='k', markeredgewidth=0,
                ms=4, zorder=0, label = 'non-proteomic data')

ax3.plot((R_P_for/2.1)*func(lambda_for, *popt_dai),lambda_for, 'o', color=  'k', #color=  '#1F4B99',
                alpha=0.2, markeredgecolor='k', markeredgewidth=0,
                ms=4, zorder=0, label = 'non-proteomic data')

ax3.plot((R_P_brem/2.1)*func(lambda_brem, *popt_dai),lambda_brem, 'o', color=  'k', #'#B7741A',
                alpha=0.2, markeredgecolor='k', markeredgewidth=0,
                ms=4, zorder=0, label = 'non-proteomic data')


for g, d in df_ribo_frac.groupby(['dataset', 'condition', 'growth_rate_hr']):
    # if g[0] == 'peebo_2015':
    #     continue
    # if g[0] == 'valgepea_2013':
    #     continue
    ax3.plot(d['frac_ribo']*func( d.growth_rate_hr.unique(), *popt_dai), g[2], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'nutrient conditions':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    ax3.plot((d['RNA/protein']/2.1)*func( d.growth_rate_hr.unique(), *popt_dai),d.growth_rate_hr.unique(), 'o', color=  'k', #color= k,
                    alpha=0.2, markeredgecolor='k', markeredgewidth=0,
                    ms=4, zorder=0, label = 'non-proteomic data') #label = g[2], )


ax3.xaxis.set_tick_params(labelsize=5)
ax3.yaxis.set_tick_params(labelsize=5)

ax3.set_xlim([0, 0.26])
ax3.set_ylim([0, 2.2])
ax3.set_xlabel('estimated active ribosomal\nfraction ($\Phi_R f_a$)', fontsize=6)
ax3.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot the prediction.
frac = np.linspace(0,0.3,100)
gr = (17.1 * frac/ L_R) * 3600
ax3.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)

axins = ax3.inset_axes([0.12, 0.65, 0.32, 0.32])
axins.plot(dai_nut_df.growth_rate_hr, dai_nut_df.f_a, 'o', color=colors['light_yellow'],
                alpha=0.9, markeredgecolor='k', markeredgewidth=0.25,
                 ms=3, zorder=1)

gr = np.linspace(0,2,100)
axins.plot(gr, func(gr, *popt_dai), color='k',
                alpha=0.5, zorder=0, lw = 0.5, ls = '--')

axins.xaxis.set_tick_params(labelsize=4)
axins.set_xticks([0,0.5,1.0,1.5,2.0])
axins.yaxis.set_tick_params(labelsize=4)
axins.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
axins.set_ylabel('$f_a$', fontsize=6)


# sub region of the original image
# axins.set_xlim(x1, x2)
# axins.set_ylim(y1, y2)
# axins.set_xticklabels('')
# axins.set_yticklabels('')

# ax3.indicate_inset_zoom(axins, alpha = 0.25, lw = 0.5)

# %
############################
# Plot 4 - rRNA limitations
############################

# Load experimental data
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
data = data[data['shorthand']=='ribosome']

# Define constants
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
n_ori = constants['N_ori']['value']

# Compute the maximum number of rRNA produced per doubling time
k_rrna = 1 # functional rRNA units per second at steady state
n_operon = 7 # average number of functional ribosomal operons per chromosome
n_rRNA_full = n_operon * k_rrna * n_ori * t_double
n_rRNA_noparallel = n_operon * k_rrna * t_double

ax4.set_yscale('log')
ax4.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax4.set_ylabel(r'number of ribosomal units', fontsize=6)
ax4.set_ylim([5E3, 8E5])
ax4.xaxis.set_tick_params(labelsize=5)
ax4.yaxis.set_tick_params(labelsize=5)

# Plot the predicted maximal number of rRNA under different regimes
ax4.plot(growth_rate, n_rRNA_full, '-', lw=1, color=colors['blue'],
        label='multiple DNA initiations\nper cell cycle (rRNA units)')
ax4.plot(growth_rate, n_rRNA_noparallel, '--', lw=1, color=colors['blue'],
        label='single DNA initiations \nper cell cycle (rRNA units)')


# Plot the experimentally observed number of ribosomes
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax4.plot(d['growth_rate_hr'], d['n_complex'], 'o',  ms=4, markeredgecolor='k',
            markeredgewidth=0.25, color=dataset_colors[g[0]], label='__nolegend__')

ax4.legend(fontsize=6, bbox_to_anchor=(1, 1))


plt.tight_layout()
plt.savefig('../../figures/fig10_ribosome_as_limit.pdf')
