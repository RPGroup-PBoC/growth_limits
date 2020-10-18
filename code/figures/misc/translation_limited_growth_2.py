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

# %%
######################
# plot 1  #
######################
# plot of avg copy num vs loc
# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# load in position information
regulonDB = pd.read_csv('../../data/regulonDB_raw/GeneProductSet.txt', delimiter = '	')

tss_map = dict(zip(regulonDB.b_number.values, regulonDB['Gene left end position in the genome'].values))
data['pos'] =data['b_number'].map(tss_map)

data_schmidt = data[data.dataset == 'schmidt_2016']
data_schmidt = data_schmidt.sort_values(by='growth_rate_hr', ascending = True)
###################
# Deal with double copies of EF-Tu
# Assume that they are equally distributed between
# gene copies. (right now, only tufA)
data_tufA = data_schmidt[data_schmidt.gene_name == 'tufA']
data_tufA['tot_per_cell'] = data_tufA['tot_per_cell']/2

data_tufB = data_schmidt[data_schmidt.gene_name == 'tufA']
data_tufB['gene_name'] = data_tufB['gene_name'].replace('tufA', 'tufB')
data_tufB['tot_per_cell'] = data_tufB['tot_per_cell']/2
data_tufB['pos'] = 4175944

data_schmidt = data_schmidt[data_schmidt.gene_name != 'tufA']
data_schmidt = data_schmidt.append(data_tufA)
data_schmidt = data_schmidt.append(data_tufB)

###################
###################

colors_viridis = plt.cm.cividis(np.linspace(0,1,len(data_schmidt.condition.unique())))
colordic = dict(zip(data_schmidt.condition.unique(), colors_viridis))


data_schmidt_ = data_schmidt[data_schmidt.pos >= (4.6E6-0.5E6)]
data_schmidt_['pos'] = data_schmidt_['pos'] - 4.6E6

data_schmidt = data_schmidt.append(data_schmidt_)

data_schmidt_ = data_schmidt[data_schmidt.pos <= (0.5E6)]
data_schmidt_['pos'] = data_schmidt_['pos'] + 4.6E6

data_schmidt = data_schmidt.append(data_schmidt_)

# an array of parameters, each of our curves depend on a specific
# value of parameters
parameters = data_schmidt.growth_rate_hr.unique()

# norm is a class which, when called, can normalize data into the
norm = matplotlib.colors.Normalize(
    vmin=np.min(parameters),
    vmax=np.max(parameters))

# choose a colormap
c_m = matplotlib.cm.cividis_r
# c_m = matplotlib.cm.magma_r

s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)

for c, d in data_schmidt.groupby('condition', sort= False):
    # pos = np.linspace(0-250000,(4.6E6-250000), 500)
    pos = np.linspace(0,(4.6E6), 500)
    avg_num = []
    for p in pos:
        # if p == 0:
        #     continue
        d_ = d[d.pos <= p + 250000]
        d_ = d_[d_.pos >= p - 250000]
        avg_num = np.append(avg_num,d_.tot_per_cell.mean())

    avg_num = avg_num - np.mean(avg_num)

    ax5.plot((pos)/1E6, avg_num,
            color = s_m.to_rgba(d.growth_rate_hr.unique())[0],
            lw = 0.5)


divider = make_axes_locatable(ax5)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(s_m, cax=cax, orientation='vertical');

cb.ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
cb.ax.tick_params(labelsize=5)
ax5.set_xlabel('genomic position (Mb)', fontsize=6)
ax5.set_ylabel('average protein copy\nnumber relative to mean', fontsize=6)
ax5.xaxis.set_tick_params(labelsize=5)
ax5.yaxis.set_tick_params(labelsize=5)
# ax5.set_xlim(0,4.6E6)
ax5.set_xlim(0,4.6)

# add in position info for rRNA, r-protein, oriC, ter
ax_genes = divider.append_axes('bottom', size='5%', pad=0.05)
ax_genes.set_xticklabels([])
ax_genes.set_yticklabels([])

ax_genes.set_xlim(0,4.6)
ax_genes.grid(False)
ax_genes.spines['bottom'].set_visible(False)
ax_genes.spines['right'].set_visible(False)
ax_genes.spines['top'].set_visible(False)
ax_genes.spines['left'].set_visible(False)
ax_genes.set_facecolor('white')
ax_genes.patch.set_alpha(0.0)

ax_genes.plot(np.array([1339769,1339769])/1E6, [-2,2], lw=2, color="#C24E9D", zorder = 10)
ax_genes.plot(np.array([1607181,1607181])/1E6, [-2,2], lw=2, color="#C24E9D", zorder = 10)
ax_genes.plot(np.array([3926090, 3926090])/1E6, [-2,2], lw=2, color="#19733A", zorder = 10)

rrna_operons = np.array([4035239, 4166367, 3941516, 3429047, 4207863, 2731448, 223593])
for p in rrna_operons:
    ax_genes.plot(np.array([p,p])/1E6, [-2,2], lw=1, color="#E55E68")

pos_schmidt = data[data.dataset == 'schmidt_2016']
data['pos'] =data['b_number'].map(tss_map)
ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']
# ribosome_genes = ['tufA']

pos_schmidt = pos_schmidt[pos_schmidt.gene_name.isin(ribosome_genes)].pos.unique()
for p in pos_schmidt:
    ax_genes.plot(np.array([p,p])/1E6, [-2,2], lw=0.5, color="#F47530")
# %%
######################
# plot 2  #
######################
# plot of t_cyc vs. tau

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


# t_cyc = data_si_mean[data_si_mean['type of perturbation'] == 'nutrient conditions']['C+D period (minutes)'].values
# t_C = data_si_mean[data_si_mean['type of perturbation'] == 'nutrient conditions']['C period (minutes)'].values
# lambda_si = data_si_mean[data_si_mean['type of perturbation'] == 'nutrient conditions']['growth_rate_hr'].values
#
# popt_tcyc, pcov_tcyc = curve_fit(func, lambda_si, t_cyc, p0=(1, 1))
# popt_tC, pcov_tC = curve_fit(func, lambda_si, t_C, p0=(1, 1))

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

# %%
######################
# plot 2  #
######################
# plot of RNA/protein vs. num ori/num ter

# Now plot!
for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'nutrient conditions':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    if 60*(np.log(2)/d.growth_rate_hr.unique()) <= 40:
        tC_dai = t_C_const
    else:
        tC_dai = func_lin(d.growth_rate_hr.unique(), *popt_tC_lin)

    ori_ter_dai = 2**(tC_dai /( np.log(2)/(d.growth_rate_hr.unique()/60)))
    ax3.plot(ori_ter_dai,  d['RNA/protein']/2.1, 'o', color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)


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


for i, l_dai in enumerate(lambda_dai):
    if 60*(np.log(2)/l_dai) <= 40:
        tC_dai = t_C_const
    else:
        tC_dai = func_lin(l_dai, *popt_tC_lin)
    # tC_dai = func(lambda_dai, *popt_tC)
    ori_ter_dai = 2**(tC_dai /( np.log(2)/(l_dai/60)))
    ax3.plot(ori_ter_dai,  R_P_dai[i]/2.1, 'o', color= colors['pale_yellow'],
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

# add in data from Scott et al.
[lambda_scott, R_P_scott] = np.array(\
        [[0.4, 0.177],
        [0.57, 0.230],
        [0.71, 0.221],
        [1.00, 0.287],
        [1.31, 0.414],
        [1.58, 0.466]]).T

for i, l_scott in enumerate(lambda_scott):
    if 60*(np.log(2)/l_scott) <= 40:
        tC_scott = t_C_const
    else:
        tC_scott = func_lin(l_scott, *popt_tC_lin)
# tC_scott = func(lambda_scott, *popt_tC)
    ori_ter_scott = 2**(tC_scott /( np.log(2)/(l_scott/60)))
    ax3.plot(ori_ter_scott,  R_P_scott[i]/2.1, 'o', color= colors['light_purple'],
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

# Forchhammer & Lindahl
[lambda_for, R_P_for] = np.array(\
        [[0.38, 0.189],
        [0.60, 0.224],
        [1.04, 0.295],
        [1.46, 0.421],
        [1.73, 0.469]]).T

for i, l_for in enumerate(lambda_for):
    if 60*(np.log(2)/l_for) <= 40:
        tC_for = t_C_const
    else:
        tC_for = func_lin(l_for, *popt_tC_lin)

    # tC_for = func(lambda_for, *popt_tC)
    ori_ter_for = 2**(tC_for /( np.log(2)/(l_for/60)))
    ax3.plot(ori_ter_for,  R_P_for[i]/2.1, 'o', color= '#1F4B99',
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

# Bremmer and Dennis
[lambda_brem, R_P_brem] = np.array(\
        [[0.42, 0.200],
        [0.69, 0.225],
        [1.04, 0.331],
        [1.39, 0.391],
        [1.73, 0.471]]).T

for i, l_brem in enumerate(lambda_brem):
    if 60*(np.log(2)/l_brem) <= 40:
        tC_brem = t_C_const
    else:
        tC_brem = func_lin(l_brem, *popt_tC_lin)
    # tC_brem = func(lambda_brem, *popt_tC)
    ori_ter_brem = 2**(tC_brem /( np.log(2)/(l_brem/60)))
    ax3.plot(ori_ter_brem,  R_P_brem[i]/2.1, 'o', color= '#B7741A',
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

ax3.set_xlabel('estimated # ori / # ter', fontsize=6)
ax3.set_ylabel('RNA/protein ratio', fontsize=6)
ax3.xaxis.set_tick_params(labelsize=5)
ax3.yaxis.set_tick_params(labelsize=5)
# ax3.set_xlim(1,4)
ax3.set_ylim(0,0.3)

# %%
######################
# plot 2  #
######################
# plot of RNA/protein vs. num ori
# Load the complex subunit counts.
subunits = pd.read_csv('../../data/compiled_annotated_complexes.csv')

# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')
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
#
# # Parameters for plotting
# rt_max = 17.1
# Kd = 15*((rt_max/ 11.5 ) - 1) #mM
# R = 12000 #copies per cell
# # ribo_frac = 0.1
# ribo_frac = 0.13
# L_R = (7459.0 + 2*394)
L_R = (7459.0)
#
# # draw colormap according to constant values of r_aa
# for raa in np.logspace(3, 9,1000):
#
#     rt = np.linspace(8,17,500)
#     yy = raa * np.ones(500)
#     fa = np.array([np.min([raa * (rt_max/rt_ - 1) / (rt_ * R * Kd), 1.0]) for rt_ in rt])
#
#     # doubling time (s):
#     tau = (L_R/(ribo_frac*fa)) * 2 * (Kd*R*fa/yy)/(np.sqrt(4*Kd*R*fa*rt_max/yy + 1) - 1)
#     # growth rate (hr-1):
#     zz = np.log(2)/(tau/3600)
#
#     # Create a set of line segments so that we can color them individually
#     points = np.array([fa*R, rt]).T.reshape(-1, 1, 2)
#     segments = np.concatenate([points[:-1], points[1:]], axis=1)
#
#     # Create a continuous norm to map from data points to colors
#     norm = plt.Normalize(0, 0.5)
#     lc = LineCollection(segments, cmap='viridis', norm=norm)
#     # Set the values used for the colormap and add line
#     lc.set_array(zz)
#     lc.set_linewidth(0.5)
#     line = ax6.add_collection(lc)
#
# for raa in np.logspace(5, 7,5):
#
#     rt = np.linspace(8,17,500)
#     yy = raa * np.ones(500)
#     fa = np.array([np.min([raa * (rt_max/rt_ - 1) / (rt_ * R * Kd), 1.0]) for rt_ in rt])
#
#     ax6.plot(fa*R, rt, zorder=10,
#               ls = '--', color = 'k', alpha = 0.5, lw = 0.5)
#
# for cond, d in dai_nut_df.groupby('condition'):
#     if d.growth_rate_hr.unique() > 0.5:
#         continue
#     ax6.plot(d['R']*d['f_a'], d['Translational elongation rate (aa/s)'], 'o',
#                 color= colors['pale_yellow'],
#                     alpha=0.8, markeredgecolor='k', markeredgewidth=0.25,
#                     ms=4, zorder=10)
#
# divider6 = make_axes_locatable(ax6_)
# cax6 = divider6.append_axes('right', size='5%', pad=0.05)
# cb6 = fig.colorbar(line, cax=cax6, orientation='vertical');
# cb6.ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
# cb6.ax.tick_params(labelsize=5)
#
# ax6.set_ylim(8,17)
# ax6.set_xlim(0,12000)
# ax6.set_ylabel('elongation rate [aa/s]', fontsize=6)
# ax6.set_xlabel('active ribosomes ($R \cdot f_a$)', fontsize=6)
# ax6.xaxis.set_tick_params(labelsize=5)
# ax6.yaxis.set_tick_params(labelsize=5)



def rt_func(K_d, raa, R, fa):
    rt_max = 17.1
    k = (6.022E23*0.6E-15*K_d)/1E9
    c = raa/(R*fa)
    return (np.sqrt(c**2 + 4*c*k*rt_max - 2*c*rt_max + rt_max**2) - c - rt_max)/(2*(k-1))


R = 13000
ribo_frac = 0.06
L_R = 7459.0
K_d = ((17.1/12.5 - 1) * 15)*1E-3

for raa in np.linspace(1E3,1E5, 50):
    raa_ = raa * np.ones(500)
    f_a = np.linspace(0.01,1,500)

    rt = rt_func(K_d, raa_, R, f_a)

    f_a, rt = np.meshgrid(f_a, rt)

    lambda_ = 3600*rt*ribo_frac*f_a/L_R

    c = ax6.contourf(R*f_a, rt, lambda_, 100, cmap = 'viridis')#,
#                    vmin=0, vmax=0.3)

for raa in np.logspace(3,6, 10):
    print(raa)
    raa_ = raa * np.ones(500)
    f_a = np.linspace(0.01,1,500)
    rt = rt_func(K_d, raa_, R, f_a)

    ax6.plot(f_a*R, rt, zorder=10,
              ls = '--', color = 'k', alpha = 0.5, lw = 0.5)


for cond, d in dai_nut_df.groupby('condition'):
    if d.growth_rate_hr.unique() > 0.5:
        continue
    ax6.plot(d['R']*d['f_a'], d['Translational elongation rate (aa/s)'],'o', color=colors['light_yellow'],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=1)

divider6 = make_axes_locatable(ax6_)
cax6 = divider6.append_axes('right', size='5%', pad=0.05)
cb6 = fig.colorbar(c, cax=cax6, orientation='vertical');
cb6.ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
cb6.ax.tick_params(labelsize=5)

ax6.set_ylim(8,17)
ax6.set_xlim(0,12000)
ax6.set_ylabel('elongation rate [aa/s]', fontsize=6)
ax6.set_xlabel('active ribosomes ($R \cdot f_a$)', fontsize=6)
ax6.xaxis.set_tick_params(labelsize=5)
ax6.yaxis.set_tick_params(labelsize=5)

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
    # ax7.plot(d['frac_active_ribo']/np.log(2), d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
    #                 alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
    #                 label = g[2], ms=4, zorder=10)
    ax7.plot(d['frac_active_ribo'], d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

for g, d in ribo_frac_act_df.groupby(['dataset', 'condition', 'growth_rate_hr']):
    if g[0] == 'peebo_2015':
        continue
    if g[0] == 'valgepea_2013':
        continue

    if 60*(np.log(2)/g[2]) <= 40:
        tC_schmidt = t_C_const
    else:
        tC_schmidt = func_lin(g[2], *popt_tC_lin)
    # tC_schmidt = func(g[2], *popt_tC)
    ori_ter_schmidt = 2**(tC_schmidt /( np.log(2)/(g[2]/60)))
    ax3.plot(ori_ter_schmidt,  d['frac_ribo'], 'o', color=dataset_colors[g[0]],
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


# plt.tight_layout()
fig.savefig('../../figures/translation_limit.pdf')
#
# ##############################
# ##############################
# ##### Additional figure for SI
#
# # add Dai et al. Cm data from their SI PDF:
# dai_Cm_df = pd.DataFrame({'condition' : ['20 mM potassium aspartate' for i in np.arange(3)],
#           'Cm (μM)' : [0, 2, 4],
#           'growth_rate_hr' : [0.33, 0.24, 0.17],
#           'Translational elongation rate (aa/s)' : [12.0 , 15.3, 15.8],
#      'f_a' : [0.756, 0.282, 0.158],
#      'RNA_P_ratio' : [0.152, 0.231, 0.292],
#     'type' : ['antibiotic' for i in np.arange(3)]},
#             columns = ['condition', 'Cm (μM)','growth_rate_hr',
#                     'Translational elongation rate (aa/s)','f_a',
#                        'RNA_P_ratio' , 'type' ])
#
#
# dai_Cm_df = dai_Cm_df.append(pd.DataFrame({'condition' : ['NQ1261(ΔptsG) in 0.2% glucose +10 mM NH4Cl' for i in np.arange(3)],
#       'Cm (μM)' : [0, 2, 3],
#       'growth_rate_hr' : [0.38, 0.16, 0.092],
#       'Translational elongation rate (aa/s)' : [12.4, 14.3, 15.2],
#  'f_a' : [0.790, 0.211, 0.099],
#  'RNA_P_ratio' : [0.160, 0.218, 0.252],
# 'type' : ['antibiotic' for i in np.arange(3)]},
#         columns = ['condition', 'Cm (μM)','growth_rate_hr',
#                 'Translational elongation rate (aa/s)','f_a',
#                    'RNA_P_ratio' , 'type' ]) )
#
#
#
# dai_Cm_df = dai_Cm_df.append(pd.DataFrame({'condition' : ['RDM+0.2% glucose +10 mM NH4Cl' for i in np.arange(3)],
#       'Cm (μM)' : [0, 4, 8],
#       'growth_rate_hr' : [1.8, 1.08, 0.57],
#       'Translational elongation rate (aa/s)' : [16.7 , 16.8, 17.3],
#  'f_a' : [0.958, 0.505, 0.243],
#  'RNA_P_ratio' : [0.476, 0.551, 0.621],
# 'type' : ['antibiotic' for i in np.arange(3)]},
#         columns = ['condition', 'Cm (μM)','growth_rate_hr',
#                 'Translational elongation rate (aa/s)','f_a',
#                    'RNA_P_ratio' , 'type' ]) )
#
# dai_Cm_df = dai_Cm_df.append(pd.DataFrame({'condition' : ['60 mM acetate +10 mM NH4Cl' for i in np.arange(3)],
#       'Cm (μM)' : [0, 3, 6],
#       'growth_rate_hr' : [0.46, 0.25, 0.18],
#       'Translational elongation rate (aa/s)' : [12.6  , 14.5, 15.6],
#  'f_a' : [0.879, 0.303, 0.170],
#  'RNA_P_ratio' : [0.172, 0.246, 0.304],
# 'type' : ['antibiotic' for i in np.arange(3)]},
#         columns = ['condition', 'Cm (μM)','growth_rate_hr',
#                 'Translational elongation rate (aa/s)','f_a',
#                    'RNA_P_ratio' , 'type' ]) )
#
# dai_Cm_df = dai_Cm_df.append(pd.DataFrame({'condition' : ['0.2% fructose +10 mM NH4Cl' for i in np.arange(3)],
#       'Cm (μM)' : [0, 4, 8],
#       'growth_rate_hr' : [0.69, 0.35, 0.21],
#       'Translational elongation rate (aa/s)' : [14.7 , 15.9 , 16.3],
#  'f_a' : [0.888, 0.290, 0.124],
#  'RNA_P_ratio' : [0.217, 0.323, 0.457],
# 'type' : ['antibiotic' for i in np.arange(3)]},
#         columns = ['condition', 'Cm (μM)','growth_rate_hr',
#                 'Translational elongation rate (aa/s)','f_a',
#                    'RNA_P_ratio' , 'type' ]) )
#
# dai_Cm_df = dai_Cm_df.append(pd.DataFrame({'condition' : ['0.2% glucose +10 mM NH4Cl' for i in np.arange(6)],
#       'Cm (μM)' : [0, 2, 4, 6, 8, 9],
#       'growth_rate_hr' : [0.98, 0.71, 0.53, 0.41, 0.33, 0.26],
#       'Translational elongation rate (aa/s)' : [15.9, 16.0, 16.1, 16.2, 16.5, 16.6],
#  'f_a' : [0.865, 0.519, 0.318, 0.222, 0.164, 0.123],
#  'RNA_P_ratio' : [0.294, 0.358, 0.440, 0.487, 0.511, 0.569],
# 'type' : ['antibiotic' for i in np.arange(6)]},
#         columns = ['condition', 'Cm (μM)','growth_rate_hr',
#                 'Translational elongation rate (aa/s)','f_a',
#                    'RNA_P_ratio' , 'type' ]) )
#
# dai_Cm_df['Naa'] = [prot.size.lambda2P(l) * 1E-15 * 6.022E23 / 110 for l in dai_Cm_df['growth_rate_hr'].values]
# dai_Cm_df['R'] = (dai_Cm_df['RNA_P_ratio']/2.1) * dai_Cm_df['Naa'] / 7459.0
#
# dai_Cm_df['ribo_frac'] = (dai_Cm_df['RNA_P_ratio']/2.1)
# dai_Cm_df = dai_Cm_df.sort_values(by = 'growth_rate_hr', ascending = False)
#
# # parameters
# rt_max = 17.1
# Kd = 15*((rt_max/ 11.5 ) - 1) #mM
# R = 120000
# raa = 1E6 # mM/s
#
# # colors, palette = prot.viz.bokeh_theme()
# color_dict = dict(zip(dai_Cm_df.condition.unique(), palette))
#
# fig2, ax_ = plt.subplots(1, 1, figsize = (8,5))
#
# for raa in [2154434,
#         5994842*2,
#         46415888/2,
#         359381366]:
#
#     rt = np.linspace(8,17,500)
#     yy = raa * np.ones(500)
#     fa = np.array([np.min([raa * (rt_max/rt_ - 1) / (rt_ * R * Kd), 1.0]) for rt_ in rt])
#
#     ax_.plot(fa*R, rt, zorder=10,
#               ls = '--', color = 'k', alpha = 0.5)
#
# for cond, d in dai_Cm_df.groupby(['condition', 'Cm (μM)'], sort = False):
#     ax_.scatter(d['R']*d['f_a'], d['Translational elongation rate (aa/s)'], zorder = 10,
#               color =  color_dict[cond[0]], label = d['condition'].unique()[0])
#
# ax_.set_ylim(8,17.5)
# ax_.set_xlim(0,119000)
# ax_.set_xlabel('active ribosomes ($R \cdot f_a$)', fontsize=16)
# ax_.set_ylabel('elongation rate [aa/s]', fontsize=16)
# ax_.xaxis.set_tick_params(labelsize=14)
# ax_.yaxis.set_tick_params(labelsize=14)
#
# handles, labels = ax_.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
#
# # Shrink current axis by 20%
# box = ax_.get_position()
# ax_.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#
# # Put a legend to the right of the current axis
# ax_.legend(by_label.values(), by_label.keys(), fontsize=12,
#     loc='center left', bbox_to_anchor=(1, 0.5))
#
#
# # plt.tight_layout()
# fig2.savefig('../../figures/SI_Dai_Cm.pdf', bbox_inches = 'tight')


#####################


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

# plot of RNA/protein vs. num ori/num ter
# additional plot of Si chlor data
fig3, ax_3 = plt.subplots(1, 1, figsize = (8,5))
# Now plot!
# for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
data_si_mean_ = data_si_mean[data_si_mean['type of perturbation'] == 'chloramphenicol']
data_si_mean_['concentration'] = data_si_mean_['concentration'].values.astype(np.float)
data_si_mean_ = data_si_mean_.sort_values(by='concentration')
alpha_ = 0.1
len = len(data_si_mean[data_si_mean['type of perturbation'] == 'chloramphenicol'].concentration.unique())
add_ = 0.9/len
count = 0
for c, d in data_si_mean_.groupby(['type of perturbation', 'concentration', 'strain'], sort=False):
    #
    if c[0] != 'chloramphenicol':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    if c[1] > 6.0:
        continue
    if c[1] <= 2.0:
        continue
    tau_C = d['C period (minutes)' ].values
    # tau = d['doubling time (minutes)'].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori_ter = 2**(tau_C /tau)
    ax_3.plot(ori_ter,  d['RNA/protein']/2.1, 'o', color= 'k',
                    alpha=alpha_+count*add_, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10, label = [c[1], c[2]])
    count +=1
    print(c[1],alpha_+count*add_)
ax_3.set_xlabel('estimated # ori / # ter', fontsize=6)
ax_3.set_ylabel('ribosomal fraction', fontsize=6)
ax_3.xaxis.set_tick_params(labelsize=5)
ax_3.yaxis.set_tick_params(labelsize=5)
# ax_3.legend()
ax_3.set_xlim(1,4)
ax_3.set_ylim(0,0.5)

fig3.savefig('../../figures/Si_Cm_data.pdf', bbox_inches = 'tight')


# plot of RNA/protein vs. num ori/num ter
# additional plot of Si chlor data
fig4, ax_4 = plt.subplots(1, 1, figsize = (8,5))
# Now plot!
for c, d in data_si_mean.groupby(['type of perturbation', 'growth media', 'strain']):
    if c[0] != 'chloramphenicol':
        continue
    if 'MG1655' in c[2]:
        k = colors['pale_red']
    elif 'NCM3722' in c[2]:
        k = colors['light_green']

    tau_cyc = d['C+D period (minutes)' ].values
    # tau = d['doubling time (minutes)'].values
    tau = 60 * (np.log(2)/d['growth_rate_hr'].values)
    ori = 2**(tau_cyc /tau)
    ax_4.plot(ori,  d['RNA/protein']/2.1, 'o', #color= k,
                    alpha=1, markeredgecolor='k', markeredgewidth=0.25,
                    ms=4, zorder=10)

ax_4.set_xlabel('estimated # ori ', fontsize=6)
ax_4.set_ylabel('RNA/protein ratio', fontsize=6)
ax_4.xaxis.set_tick_params(labelsize=5)
ax_4.yaxis.set_tick_params(labelsize=5)
# ax_4.set_xlim(1,4)
# ax_4.set_ylim(0,0.3)

fig4.savefig('../../figures/Si_Cm_data_2.pdf', bbox_inches = 'tight')
