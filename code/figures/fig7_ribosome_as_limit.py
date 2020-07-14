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

# Add in data from Dai et al.
[lambda_dai, R_P_dai] = np.array(\
        [[0.0, 0.08853279947047754],
        [0.03324706890845652, 0.09742834356027935],
        [0.12844176066233703, 0.12157153165470358],
        [0.19652012565308674, 0.12917148174069676],
        [0.23055930814846148, 0.13297145678369332],
        [0.2849547247405284, 0.14694954887503597],
        [0.33601892391911736, 0.15201256530868013],
        [0.4074068045812377, 0.17107537557577435],
        [0.417639175984852, 0.16979497279144085],
        [0.4517000602223341, 0.17104716331103476],
        [0.485674137491387, 0.18249049192423916],
        [0.5503561798423366, 0.1888187199227418],
        [0.6727865579409387, 0.21549233114688282],
        [0.6864152519843529, 0.21548365045003987],
        [0.7000547968988209, 0.21420107749149564],
        [0.7170798135820351, 0.21546411888214323],
        [0.744196140345166, 0.2320073568905744],
        [0.9177883754618401, 0.25227895419304786],
        [0.9448830004828637, 0.27136997672488156],
        [0.9926268331190284, 0.2662440252391261],
        [0.9753630972726335, 0.29300661360590724],
        [1.0979236858238794, 0.3043935176896325],
        [1.1624538159800777, 0.32855623735195344],
        [1.2677832212980895, 0.36288405301735593],
        [1.566952587118931, 0.4404005056505911],
        [1.7949076862145108, 0.4784718718295111]]).T

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

# inset axes....
# axins = ax2.inset_axes([0.5, 0.5, 0.47, 0.47])
# axins = ax2.inset_axes([0.5, 0.05, 0.5, 0.5])

# axins.imshow(Z2, extent=extent, interpolation="nearest",
#           origin="lower")
#
# for ax_ in [ax2, axins]:
#
#     for g, d in df_ribo_frac.groupby(['dataset', 'condition', 'growth_rate_hr']):
#         # if g[0] == 'peebo_2015':
#         #     continue
#         # if g[0] == 'valgepea_2013':
#         #     continue
#         ax_.plot(d['frac_ribo'], g[2], 'o', color=dataset_colors[g[0]],
#                         alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
#                         label = g[2], ms=4, zorder=10)
#
#
#     ax_.plot(R_P_dai/2.1, lambda_dai, 'o', color= colors['pale_yellow'],
#                         alpha=1, markeredgecolor='k', markeredgewidth=0.25,
#                         ms=4, zorder=10)
#
#     ax_.plot(R_P_scott/2.1, lambda_scott, 'o', color= colors['light_purple'],
#                     alpha=1, markeredgecolor='k', markeredgewidth=0.25,
#                     ms=4, zorder=10)
#
#
#     ax_.plot(R_P_for/2.1, lambda_for, 'o', color= '#1F4B99',
#                     alpha=1, markeredgecolor='k', markeredgewidth=0.25,
#                     ms=4, zorder=10)
#
#
#     ax_.plot(R_P_brem/2.1, lambda_brem, 'o', color= '#B7741A',
#                     alpha=1, markeredgecolor='k', markeredgewidth=0.25,
#                     ms=4, zorder=10)
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

## grab nutrient limitation data from Dai et al. 2016
# nutrient limitation data from SI PDF:
conditions = ['RDM + 0.2% glucose+10 mM NH4Cl',
            '0. 2 % glucose+cAA+10 mM NH4Cl',
           ' 10 mM glucose-6-phosphate+10 mM gluconate +10 mM NH4Cl',
            '0.2% glucose+10 mM NH4Cl',
            '0.2% xylose+10 mM NH4Cl',
            '0.2 % glycerol+10 mM NH4Cl',
            '0.2% fructose+10 mM NH4Cl',
            '0.2% sorbitol+10 mM NH4Cl',
            '0.2% galactose+10 mM NH4Cl',
            '60 mM acetate+10 mM NH4Cl',
            '0.2% mannose+10 mM NH4Cl',
            '0.1% mannose+10 mM NH4Cl',
            '20 mM potassium aspartate',
            '0.075% mannose+10 mM NH4Cl',
            '20 mM aspartate+10 mM NH4Cl',
            '0.2% glycerol +10 mM Arginine',
            '20 mM glutamate+10 mM NH4Cl',
            '0.2% glycerol+20 mM Threonine']

erate = [16.7, 16.3, 16.1, 15.9, 14.9, 15.0, 14.7, 13.7, 13.1, 12.6, 13.0, 12.4,
        12.0, 12.1, 12.3, 11.6, 10.7, 9.4]

RNA_P = [0.476, 0.364, 0.306, 0.294, 0.233, 0.227, 0.217, 0.193, 0.184, 0.172,
        0.172, 0.152, 0.152, 0.147, 0.137, 0.130, 0.118, 0.097]

r_prot_frac = [np.nan, np.nan, np.nan, 11.6, np.nan, np.nan, np.nan, np.nan,
        np.nan, 7.2, np.nan, np.nan, np.nan, 5.1, np.nan, 6.1, 4.7, 4.4]


gr_dai, fa = np.array([[1.8, 0.958],
            [1.28, 0.9],
            [1.12, 0.927],
            [0.98, 0.865],
            [0.75, 0.902],
            [0.69, 0.849],
            [0.69, 0.888],
            [0.55, 0.879],
            [0.5, 0.860],
            [0.46, 0.879],
            [0.41, 0.756],
            [0.34, 0.751],
            [0.33, 0.756],
            [0.29, 0.683],
            [0.23, 0.590],
            [0.201, 0.554],
            [0.13, 0.441],
            [0.035, 0.168]]).T

dai_nut_df = pd.DataFrame({'condition' : conditions,
     'Cm (μM)' : np.zeros(len(gr_dai)),
     'RNA_P_ratio' : RNA_P,
     'growth_rate_hr' : gr_dai,
     'Translational elongation rate (aa/s)' : erate,
     'measured_prot_frac' : r_prot_frac,
     'f_a' : fa,
    'type' : ['nutrient limitation' for i in np.arange(len(gr_dai))]},
        columns = ['condition', 'Cm (μM)', 'RNA_P_ratio', 'growth_rate_hr',
                   'Translational elongation rate (aa/s)', 'measured_prot_frac',
                   'f_a', 'type'])

# fit measurements of active fraction from Dai et al. data
dai_nut_df = dai_nut_df.sort_values(by='growth_rate_hr', ascending = True)

popt_dai, pcov_dai = curve_fit(func, dai_nut_df.growth_rate_hr.values, dai_nut_df.f_a.values, p0=(1, 1e-6, 1))


# Add Dai et al, Scott, and historical
ax3.plot((R_P_dai/2.1)*func(lambda_dai, *popt_dai),lambda_dai, 'o', color=colors['light_yellow'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax3.plot((R_P_scott/2.1)*func(lambda_scott, *popt_dai),lambda_scott, 'o', color=colors['light_purple'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax3.plot((R_P_for/2.1)*func(lambda_for, *popt_dai),lambda_for, 'o', color=  '#1F4B99',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax3.plot((R_P_brem/2.1)*func(lambda_brem, *popt_dai),lambda_brem, 'o', color=  '#B7741A',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)


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

    ax3.plot((d['RNA/protein']/2.1)*func( d.growth_rate_hr.unique(), *popt_dai),d.growth_rate_hr.unique(), 'o', color= k,
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=0)


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
axins.plot(gr_dai, fa, 'o', color=colors['light_yellow'],
                alpha=0.9, markeredgecolor='k', markeredgewidth=0.25,
                 ms=3, zorder=1)

gr = np.linspace(0,2,100)
axins.plot(gr, func(gr, *popt_dai), color='k',
                alpha=0.5, zorder=0, lw = 0.5, ls = '--')

axins.xaxis.set_tick_params(labelsize=4)
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
plt.savefig('../../figures/fig7_ribosome_as_limit.pdf')
