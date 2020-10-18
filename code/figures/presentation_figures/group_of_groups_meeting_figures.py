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
def func(x, a, c, d):
    return a*np.exp(-c*x)+d
L_R = 7459.0 # aa
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
fig = plt.figure(constrained_layout=True)
# widths = [6, 2.5, 2.5, 5]
widths = [6, 5, 5]
heights = [1.5, 1, 0.75, 1.25, 0.75, 2]
spec = fig.add_gridspec(ncols=3, nrows=6, width_ratios=widths,
                          height_ratios=heights)

# plot of t_cyc vs. tau
ax1 = fig.add_subplot(spec[0, 1])
# plot of t_cyc vs. tau
ax2 = fig.add_subplot(spec[0, 2])
# plot of RNA/protein vs. num ori
ax3 = fig.add_subplot(spec[1:4, 1])
# plot of # ribosomes vs. num ori
ax4 = fig.add_subplot(spec[1:4, 2])
# plot of avg copy num vs loc
ax5 = fig.add_subplot(spec[3, 0])


# plot of elongation rate vs. active fraction
ax6 = fig.add_subplot(spec[4:, 1])
ax6_ = fig.add_subplot(spec[4:, 0])
ax6_.axis('off')

# plot of ribosomal fraction vs lambda
ax7 = fig.add_subplot(spec[4:, 2])

######################
# plot 2  #
######################
# plot of t_cyc vs. tau

# Load the full Si 2017 SI
data_si = pd.read_csv('../../../data/si_2017_raw/si2017_full.csv')

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

# Now plot!
for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'nutrient conditions':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    ax1.plot(d['growth_rate_hr'],  d['C period (minutes)'], 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)
    ax2.plot(d['growth_rate_hr'],  d['C+D period (minutes)'], 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)


tau = np.linspace(0.1, 40, 100)
x = (np.log(2)/tau)*60
tC_x = t_C_const.mean()*np.ones(100)
ax1.plot(x,  tC_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')

tau = np.linspace(40, (np.log(2)/0.1)*60, 100)
x = (np.log(2)/tau)*60
tC_x = func_lin(x, *popt_tC_lin)
ax1.plot(x,  tC_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')

# plot 2 ; t_cyc vs. tau and lambda
tau = np.linspace(0.1, 40, 100)
x = (np.log(2)/tau)*60
tcyc_x = t_cyc_const.mean()*np.ones(100)
ax2.plot(x,  tcyc_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')

tau = np.linspace(40, (np.log(2)/0.1)*60, 100)
x = (np.log(2)/tau)*60
tcyc_x = func_lin(x, *popt_tcyc_lin)
ax2.plot(x,  tcyc_x, color = 'k', lw = 0.5,
                alpha=0.75,
                zorder=1, ls = '--')
# x = np.linspace(0.1,2,100)
# tcyc_x = func(x, *popt_tcyc)
# ax2.plot(x,  tcyc_x, color = 'k', lw = 0.5,
#                 alpha=0.75,
#                 zorder=1, ls = '--')

for a in [ax1,ax2]:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
    a.set_xlim([0, 2])
    a.set_ylim([0, 150])

ax2.set_ylabel('t$_{cyc}$ [min]', fontsize=6)
ax1.set_ylabel('t$_{C}$ [min]', fontsize=6)


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

# Forchhammer & Lindahl
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



# %%
######################
# plot 2  #
######################
# plot of RNA/protein vs. num ori
# Load the complex subunit counts.
subunits = pd.read_csv('../../../data/compiled_annotated_complexes.csv')

# # Load the compiled data
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
# data['gene_name'] = data['gene_name'].str.lower()

# Compute the minimum number of complexes.
complex_count = subunits.groupby(['dataset', 'dataset_name', 'condition',
                    'growth_rate_hr', 'complex_annotation',
                    'complex'])['n_units'].mean().reset_index()

complex_ribo = complex_count[complex_count.complex_annotation == 'ribosome']

# complex_ribo['t_cyc'] = func(complex_ribo['growth_rate_hr'], *popt_tcyc)

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
# print(t_cyc_const, complex_ribo['t_cyc'], complex_ribo['# ori'], complex_ribo['growth_rate_hr'])


for g, d in complex_ribo.groupby(['dataset', 'condition', 'growth_rate_hr']):
    ax4.plot(d['# ori'], d['n_units'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

ax4.set_xlabel('estimated # ori', fontsize=6)
ax4.set_ylabel('ribosomes per cell', fontsize=6)
ax4.xaxis.set_tick_params(labelsize=5)
ax4.yaxis.set_tick_params(labelsize=5)
# ax4.set_xlim(1,3.6)

# %%
######################
# plot 6  #
######################
# plot of elongation rate versus active ribosomes

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


gr, fa = np.array([[1.8, 0.958],
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
     'Cm (μM)' : np.zeros(len(gr)),
     'RNA_P_ratio' : RNA_P,
     'growth_rate_hr' : gr,
     'Translational elongation rate (aa/s)' : erate,
     'measured_prot_frac' : r_prot_frac,
     'f_a' : fa,
    'type' : ['nutrient limitation' for i in np.arange(len(gr))]},
        columns = ['condition', 'Cm (μM)', 'RNA_P_ratio', 'growth_rate_hr',
                   'Translational elongation rate (aa/s)', 'measured_prot_frac',
                   'f_a', 'type'])


dai_nut_df['Naa'] = [prot.size.lambda2P(l) * 1E-15 * 6.022E23 / 110 for l in dai_nut_df['growth_rate_hr'].values]
dai_nut_df['R'] = (dai_nut_df['RNA_P_ratio']/2.1) * dai_nut_df['Naa'] / 7459.0
dai_nut_df['R2'] = (dai_nut_df['measured_prot_frac']/100) * dai_nut_df['Naa'] / 7459.0
# ################

# %%
######################
# plot 7  #
######################
# plot of ribosomal fraction vs lambda
dai_nut_df = dai_nut_df.sort_values(by='growth_rate_hr', ascending = True)

popt_dai, pcov_dai = curve_fit(func, dai_nut_df.growth_rate_hr.values, dai_nut_df.f_a.values, p0=(1, 1e-6, 1))


ribo_frac_act_df = pd.DataFrame()
for c, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    mass_ribo = d[d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    mass_nonribo = d[~d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    frac_ribo = (mass_ribo )/ d.fg_per_cell.sum()

    lambda_max = (np.log(2) * 17.1 * frac_ribo/ L_R) * 3600

    data_list = {'frac_ribo' : frac_ribo,
                 'frac_active_ribo' : frac_ribo * np.min([func( c[2], *popt_dai),1.0]),
                'dataset' : c[0],
                'dataset_name' : d.dataset_name.values[0],
                'condition' : c[1],
                'growth_rate_hr' : c[2],
                'lambda_max' : lambda_max}

    ribo_frac_act_df = ribo_frac_act_df.append(data_list,
                                        ignore_index = True)

for g, d in ribo_frac_act_df.groupby(['dataset', 'condition', 'growth_rate_hr']):
    if g[0] == 'peebo_2015':
        continue
    if g[0] == 'valgepea_2013':
        continue

    ax7.plot(d['frac_active_ribo'], d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)




# Add Dai et al, Scott, and historical
ax7.plot((R_P_dai/2.1)*func(lambda_dai, *popt_dai),lambda_dai, 'o', color=colors['light_yellow'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                label = g[2], ms=4, zorder=1)

ax7.plot((R_P_scott/2.1)*func(lambda_scott, *popt_dai),lambda_scott, 'o', color=colors['light_purple'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                label = g[2], ms=4, zorder=1)

ax7.plot((R_P_for/2.1)*func( lambda_for, *popt_dai),lambda_for, 'o', color=  '#1F4B99',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                label = g[2], ms=4, zorder=1)

ax7.plot((R_P_brem/2.1)*func( lambda_brem, *popt_dai),lambda_brem, 'o', color=  '#B7741A',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                label = g[2], ms=4, zorder=1)

for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'nutrient conditions':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']
    ax7.plot((d['RNA/protein']/2.1)*func( d.growth_rate_hr.unique(), *popt_dai),d.growth_rate_hr.unique(), 'o', color= k,
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=0)


# add translation limited growth rate
frac = np.linspace(0,0.3,100)
# gr = (np.log(2) * 17.1 * frac/ L_R) * 3600
gr = (17.1 * frac/ L_R) * 3600
ax7.plot(frac, gr,
        lw = 0.5, color = 'k')

ax7.set_xlim(0,0.3)
ax7.set_ylim(np.min(gr), np.max(gr))



ax7.set_ylabel('growth rate [hr$^{-1}]$', fontsize=6)
ax7.set_xlabel('estimated active ribosomal\nfraction ($\Phi_R \cdot f_a$)', fontsize=6)
ax7.xaxis.set_tick_params(labelsize=5)
ax7.yaxis.set_tick_params(labelsize=5)




dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}

# # Load the compiled data
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')

# calculate ribosomal fraction
L_R = 7459.0 # aa
ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']

ribo_frac_df = pd.DataFrame()
for c, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    mass_ribo = d[d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    frac_ribo = mass_ribo / d.fg_per_cell.sum()
    data_list = {'frac_ribo' : frac_ribo,
                'dataset' : c[0],
                'condition' : c[1],
                'growth_rate_hr' : c[2],
                'dataset_name' : d.dataset_name.values[0]}
    ribo_frac_df = ribo_frac_df.append(data_list,
                                        ignore_index = True)

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 0.3])
ax.set_ylim([0, 2.2])
ax.set_xlabel('ribosomal fraction ($\Phi_R$)', fontsize=6)
ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot the prediction.
frac = np.linspace(0,1.0,100)
gr = (17.1 * frac/ L_R) * 3600
ax.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)
ax.hlines((17.1 / L_R) * 3600, 0, np.max(gr), color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in ribo_frac_df.groupby(['dataset', 'dataset_name']):
    ax.plot(d['frac_ribo'], d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[1], ms=4)

plt.savefig('../../../figures/g-of-g_ribosome_frac_vs_rate_a.pdf', bbox_inches='tight')
# Add Dai et al, Scott, and historical
ax.plot((R_P_dai/2.1),lambda_dai, 'o', color=colors['light_yellow'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_scott/2.1),lambda_scott, 'o', color=colors['light_purple'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_for/2.1),lambda_for, 'o', color=  '#1F4B99',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_brem/2.1),lambda_brem, 'o', color=  '#B7741A',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

# ax.legend(ncol=1, fontsize=6)
plt.savefig('../../figures/g-of-g_ribosome_frac_vs_rate_b.pdf', bbox_inches='tight')


#############################################
#############################################
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 0.3])
ax.set_ylim([0, 2.2])
ax.set_xlabel('ribosomal fraction ($\Phi_R$)', fontsize=6)
ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot the prediction.
frac = np.linspace(0,1.0,100)
gr = (17.1 * frac/ L_R) * 3600
ax.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)
ax.hlines((17.1 / L_R) * 3600, 0, np.max(gr), color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in ribo_frac_df.groupby(['dataset', 'dataset_name']):
    if g[0] == 'peebo_2015':
        continue
    if g[0] == 'valgepea_2013':
        continue
    ax.plot(d['frac_ribo'], d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[1], ms=4)

# Add Dai et al, Scott, and historical
ax.plot((R_P_dai/2.1),lambda_dai, 'o', color=colors['light_yellow'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_scott/2.1),lambda_scott, 'o', color=colors['light_purple'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_for/2.1),lambda_for, 'o', color=  '#1F4B99',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_brem/2.1),lambda_brem, 'o', color=  '#B7741A',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

# ax.legend(ncol=1, fontsize=6)
plt.savefig('../../../figures/g-of-g_ribosome_frac_vs_rate_b_2.pdf', bbox_inches='tight')

#############################################
#############################################
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 0.3])
ax.set_ylim([0, 2.2])
ax.set_xlabel('estimated active ribosomal\nfraction ($\Phi_R \cdot f_a$)', fontsize=6)
ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot the prediction.
frac = np.linspace(0,1.0,100)
gr = (17.1 * frac/ L_R) * 3600
ax.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)
ax.hlines((17.1 / L_R) * 3600, 0, np.max(gr), color='k', linestyle='--', lw=0.75, label='__nolegend__')


for g, d in ribo_frac_act_df.groupby(['dataset', 'condition', 'growth_rate_hr']):
    if g[0] == 'peebo_2015':
        continue
    if g[0] == 'valgepea_2013':
        continue
    ax.plot(d['frac_active_ribo'], d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[1], ms=4)

# Add Dai et al, Scott, and historical
ax.plot((R_P_dai/2.1)*func( lambda_dai, *popt_dai),lambda_dai, 'o', color=colors['light_yellow'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_scott/2.1)*func( lambda_scott, *popt_dai),lambda_scott, 'o', color=colors['light_purple'],
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_for/2.1)*func( lambda_for, *popt_dai),lambda_for, 'o', color=  '#1F4B99',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

ax.plot((R_P_brem/2.1)*func( lambda_brem, *popt_dai),lambda_brem, 'o', color=  '#B7741A',
                alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                ms=4, zorder=1)

# ax.legend(ncol=1, fontsize=6)
plt.savefig('../../../figures/g-of-g_act_ribosome_frac_vs_rate.pdf', bbox_inches='tight')


# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# # Load the compiled data
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')


data['tau']  = 60*(np.log(2)/(data['growth_rate_hr']))
t_cyc_arr = []
for i, val in enumerate(data['tau'].values):
    if val <= 40:
        t_cyc_ = t_cyc_const
    else:
        t_cyc_ = func_lin(data['growth_rate_hr'].values[i], *popt_tcyc_lin)
    t_cyc_arr = np.append(t_cyc_arr, t_cyc_)

data['t_cyc'] = t_cyc_arr #func_lin(complex_ribo['growth_rate_hr'], *popt_tcyc_lin)
data['# ori'] =  2**(data['t_cyc'] / data['tau'] )
# print(t_cyc_const, complex_ribo['t_cyc'], complex_ribo['# ori'], complex_ribo['growth_rate_hr'])


for g, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    ax.plot(d['# ori'].unique(), d['fg_per_cell'].sum(), 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

ax.set_xlabel('estimated # ori', fontsize=6)
ax.set_ylabel('total protein mass\n per cell', fontsize=6)
ax.xaxis.set_tick_params(labelsize=5)
ax.yaxis.set_tick_params(labelsize=5)

plt.savefig('../../../figures/g-of-g_ribosome_ori_vs_fg.pdf', bbox_inches='tight')






# %%
######################
# plot 2  #
######################
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# plot of RNA/protein vs. num ori
# Load the complex subunit counts.
subunits = pd.read_csv('../../../data/compiled_annotated_complexes.csv')

# # Load the compiled data
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
# data['gene_name'] = data['gene_name'].str.lower()

# Compute the minimum number of complexes.
complex_count = subunits.groupby(['dataset', 'dataset_name', 'condition',
                    'growth_rate_hr', 'complex_annotation',
                    'complex'])['n_units'].mean().reset_index()

complex_ribo = complex_count[complex_count.complex_annotation == 'ribosome']

# complex_ribo['t_cyc'] = func(complex_ribo['growth_rate_hr'], *popt_tcyc)

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
# print(t_cyc_const, complex_ribo['t_cyc'], complex_ribo['# ori'], complex_ribo['growth_rate_hr'])


for g, d in complex_ribo.groupby(['dataset', 'condition', 'growth_rate_hr']):
    ax.plot(d['# ori'], d['n_units'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

ax.set_xlabel('estimated # ori', fontsize=6)
ax.set_ylabel('ribosomes per cell', fontsize=6)
ax.xaxis.set_tick_params(labelsize=5)
ax.yaxis.set_tick_params(labelsize=5)
plt.savefig('../../../figures/g-of-g_ribosome_ori_vs_ribosomes.pdf', bbox_inches='tight')




# %%
######################
######################
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

aa = np.linspace(0.01,75)
er = 17.1/(1 + 6/aa)

ax.plot(aa,er, lw = 1, color = 'k')
ax.set_xlabel('[AA]', fontsize=6)
ax.set_ylabel('elongation rate [aa/s]', fontsize=6)
ax.xaxis.set_tick_params(labelsize=5)
ax.yaxis.set_tick_params(labelsize=5)
ax.set_xticks([])

plt.savefig('../../../figures/g-of-g_elongation_rate.pdf', bbox_inches='tight')
