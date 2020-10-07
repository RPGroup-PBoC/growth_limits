#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import prot.viz
import prot.size as size
colors, palette = prot.viz.bokeh_theme()
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
prot.viz.plotting_style()


######################
# plot configuration #
######################
fig = plt.figure(constrained_layout=True)
widths = [6, 7]
heights = [2, 0.5, 0.5, 0.25]
spec = fig.add_gridspec(ncols=2, nrows=4, width_ratios=widths,
                          height_ratios=heights)

ax1 = fig.add_subplot(spec[0, 0])
# ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[2:, 0])
ax4 = fig.add_subplot(spec[:2, 1])


# Parameters and calculations #
###############################
# use estimates of E. coli rod cell size based
# on the values used from Si et al. (2017,2019);
# Use our fit of that data as fcn of growth rate
# to determine cell size, length, width, and surface area.
# gr = np.linspace(0.01, 4.4, 100)
V = np.linspace(0.5,50, 500)
# V = size.lambda2size(gr)
# print(V)
# w = size.lambda2width(gr)
# l = size.lambda2length(gr)
# SA_rod = size.rod_SA(l, w, V)

SA_rod = 2 * np.pi *  V**(2/3)
SA_V_ratio_rod = SA_rod / V

# Sphere V = (4/3) pi r**3
# Sphere SA = 4 pi r**2
SA_sphere = V**(2/3) * ((4/3) * np.pi)**(-2/3) * 4 * np.pi
SA_V_ratio_sphere = SA_sphere / V

# ATP equivalents demand w.r.t. volume ; 1E6 ATP/(um3 s)
Pv = 1E6 * V

# ATP max - half surface area devoted to respiration
Ps_um_resp = ((3)/ (1E-6))
Ps_resp_rod = Ps_um_resp * SA_rod * 0.5
Ps_resp_sphere = Ps_um_resp * SA_sphere * 0.5

# # ATP max - half surface area devoted to fermentation
# Ps_um_ferm= ((180*2)/ (50E-6))
# Ps_ferm_rod = Ps_um_ferm * SA_rod * 0.5
# Ps_ferm_sphere = Ps_um_ferm * SA_sphere * 0.5

#### for Fill_between, need common x-axis array
SA_ = np.linspace(np.min(SA_sphere), np.max(SA_rod), 500)
V_sphere = (SA_/(((4/3) * np.pi)**(-2/3) * 4 * np.pi))**(3/2)
SA_V_ratio_sphere_ = SA_/V_sphere

V_rod = (SA_/(2 * np.pi))**(3/2)
SA_V_ratio_rod_ = SA_/V_rod

# ATP max - half surface area devoted to respiration
Ps_um_resp = ((3)/ (1E-6))
Ps_resp_ = Ps_um_resp * SA_ * 0.5



######################
# Plot 1, S/V scaling #
######################

ax1.plot(Pv, SA_V_ratio_rod, color=colors['dark_green'], label='rod',
           alpha=0.9, lw = 0.5, ls = '-.')
ax1.plot(Pv, SA_V_ratio_sphere, color=colors['dark_green'], label='sphere',
           alpha=0.9, lw = 0.5, ls = '--')

ax1.fill_between(Pv, y1 = SA_V_ratio_sphere, y2=SA_V_ratio_rod,
        color=colors['dark_green'],alpha=0.2, lw = 0)

# plot the max for respiration
ax1.plot(Ps_resp_rod, SA_V_ratio_rod, color=colors['blue'],
            label='rod', alpha=0.9, lw = 0.5, ls = '-.')
ax1.plot(Ps_resp_sphere, SA_V_ratio_sphere, color=colors['blue'],
            label='sphere', alpha=0.9, lw = 0.5, ls = '--')

ax1.fill_between(Ps_resp_, y1 = SA_V_ratio_sphere_, y2 = SA_V_ratio_rod_,
        color=colors['blue'],alpha=0.2, lw = 0)

# # Populate second plot with growth rates
# S/V for E. coli datasets
# Load the data set
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

for g, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    V = size.lambda2size(g[2])
    # ATP equivalents demand w.r.t. volume ; 1E6 ATP/(um3 s)
    Pv = 1E6 * V
    # assume aspect ratio of 4 (length/width), which is
    # appoximately correct for E. coli
    SA_rod = 2 * np.pi *  V**(2/3)
    SV = SA_rod/V

    ax1.plot(Pv, SV, 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

# Format the axes
for a in [ax1]:#,ax2]:
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
    a.set_xscale('log')
    a.set_ylim([1.5, 8.0])
    # a.legend(fontsize=5, loc='lower right')

ax1.set_xlim([np.min(Pv), np.max(Ps_resp_)])
ax1.set_xlabel('ATP equivalents per s', fontsize=6)
ax1.set_ylabel('S/V ratio [$\mu$m]', fontsize=6)



legend_elements = [Line2D([0], [0], marker='o', color='gray', lw = 0,
                   alpha=0.5, markeredgecolor='k', markeredgewidth=0.25,
                   label = 'estimated for\n$E. coli$ data', ms=4),
                   Line2D([0], [0], ls = '-.', color='gray', label='rod (aspect\nratio = 4)',
                     alpha = 1.0, lw = 0.5, markeredgewidth=0.0, markeredgecolor='k'),
                   Line2D([0], [0], ls = '--', color='gray', label='sphere',
                     alpha = 1.0, lw = 0.5, markeredgewidth=0.0, markeredgecolor='k')]

legend = ax1.legend(handles = legend_elements, loc = 'upper right',
                        fontsize = 5)
ax1.add_artist(legend)



######################
# Plot 2
######################
# total fg per cell of inner membrane proteins, GO:0005886
data_membrane = data[data.go_terms.str.contains('GO:0005886')]

data_membrane_fg_summary = pd.DataFrame()
for c, d in data_membrane.groupby(['dataset', 'condition', 'growth_rate_hr']):
    d_ = d[d.gene_name != 'tufA']
    d_ = d_[d_.gene_name != 'tufB']

    # V = 0.28 * np.exp(1.33 * c[2])
    # SA = 2 * np.pi *  V**(2/3)
    w = size.lambda2width(c[2])
    l = size.lambda2length(c[2])
    V = size.lambda2size(c[2])
    SA = size.rod_SA(l, w, V)

    fg_tot = d_.fg_per_cell.sum()
    fg_SA = fg_tot / SA

    data_list = {'fg per um2' : fg_SA,
                'growth_rate_hr' : c[2],
                'dataset' : c[0],
                'condition' : c[1]}
    data_membrane_fg_summary = \
            data_membrane_fg_summary.append(data_list,
            ignore_index = True)

for g, d in data_membrane_fg_summary.groupby(['dataset', 'condition', 'growth_rate_hr']):
    ax3.plot(g[2], d['fg per um2'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[0], ms=4, zorder=10)
ax3.set_ylim(0,10)
ax3.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax3.set_ylabel('[fg per $\mu m^2$]', fontsize=6)
ax3.xaxis.set_tick_params(labelsize=5)
ax3.yaxis.set_tick_params(labelsize=5)



######################
# Plot 3
######################
# Plot distribution of inner membrane proteins, GO:0005886
# Load the complex subunit counts.
subunits = pd.read_csv('../../data/compiled_annotated_complexes.csv')
complex_energy = ['NADH-DHI-CPLX', 'CPLX0-8160', 'CYT-O-UBIOX-CPLX', 'CYT-D-UBIOX-CPLX', 'ATPSYN-CPLX']

df_mem = data_membrane[data_membrane.dataset == 'schmidt_2016']
df_mem = df_mem[df_mem.gene_name != 'tufA']

genes_respiration = []
for c, d in subunits.groupby('complex'):
    if c in complex_energy:
        for gene, d_ in d.groupby('gene_name'):
            genes_respiration = np.append(genes_respiration, gene)

genes_carbon = subunits[subunits.go_terms.str.contains('GO:0008643')].gene_name.unique()

df_mem = df_mem.replace(genes_respiration, 'respiration')
df_mem = df_mem.replace(genes_carbon, 'carbon_uptake')

df_mem = df_mem.replace('Not Assigned', 'poorly characterized or not assigned')
df_mem = df_mem.replace('poorly characterized', 'poorly characterized or not assigned')

df_mem_ = pd.DataFrame()
for c,d in df_mem.groupby(['dataset', 'condition',
                                   'growth_rate_hr', 'gene_name']):
    data_list = {'gene_name' : d.gene_name.unique()[0] ,
                 'cog_class' : d.cog_class.unique()[0] ,
                 'condition' : d.condition.unique()[0] ,
                 'dataset' : d.dataset.unique()[0] ,
                 'dataset_name' : d.dataset_name.unique()[0] ,
                 'fg_per_cell' : d.fg_per_cell.sum(),
                'growth_rate_hr' :  d.growth_rate_hr.unique()[0]}
    df_mem_ = df_mem_.append(data_list,
                                            ignore_index=True)

df_mem_['rel_fg_per_cell'] = df_mem_.groupby('condition').transform(lambda x: (x / x.sum()))['fg_per_cell']
df_mem_ = df_mem_.sort_values(by=['growth_rate_hr', 'gene_name'], ascending = False)

y_order = dict(zip(df_mem_.condition.unique(), np.arange(len(df_mem_.condition.unique()))))
cog_class_order = ['metabolism',
                'cellular processes and signaling',
                'information storage and processing',
                'poorly characterized or not assigned']
order_dict = dict(zip(cog_class_order,
                np.arange(4)))
# color_dict = dict(zip(cog_class_order,
#                     ['', '#C8715B', '#A587AA', '#788FBD']))
color_dict = dict(zip(cog_class_order,
                    ['', '#BF703A', '#D3B15E', '#788FBD']))


for c, d in df_mem_.groupby('condition', sort=False):
    for c_ in cog_class_order:
        if c_ == 'metabolism':
            resp = d[d.gene_name == 'respiration']
            ax4.barh(y_order[c], resp.rel_fg_per_cell, height=0.9, color='#679B48', alpha=0.3, #84A779
                        linewidth=0.1)

            c_uptake = d[d.gene_name == 'carbon_uptake']
            lefts = resp.rel_fg_per_cell.sum()
            ax4.barh(y_order[c], c_uptake.rel_fg_per_cell.sum(),  height=0.9, color='#679B48', alpha=0.7,
                            left=lefts, linewidth=0.1)

            # lefts = d[d.gene_name == 'respiration']
            lefts += d[d.gene_name == 'carbon_uptake'].rel_fg_per_cell.sum()
            meta = d[d.gene_name != 'respiration'].copy()
            meta = meta[meta.gene_name !='carbon_uptake'].copy()

            ax4.barh(y_order[c], meta.rel_fg_per_cell.sum(), height=0.9, color='#679B48',
                            left=lefts,  linewidth=0.1)

            lefts += meta[meta.cog_class == c_].rel_fg_per_cell.sum()

        else:
            ax4.barh(y_order[c], d[d.cog_class == c_].rel_fg_per_cell.sum(), height=0.9,
                        color=color_dict[c_], left=lefts, linewidth=0.1)
            lefts += d[d.cog_class == c_].rel_fg_per_cell.sum()



ax4.set_xlim(0,1)
ax4.set_ylim(0,len(df_mem_.condition.unique()))
ax4.set_yticks(np.arange(len(df_mem_.condition.unique()))-0.5)
ax4.set_yticklabels(df_mem_.condition.unique())
ax4.set_xlabel('relative plasma membrane abundance (GO term : 0005886)', fontsize=6)
ax4.xaxis.set_tick_params(labelsize=5)
ax4.yaxis.set_tick_params(labelsize=5)

growth_rates_list = [d.growth_rate_hr.unique()[0] for c, d in df_mem_.groupby('condition', sort=False)]
ax4_twin = ax4.twinx()
ax4_twin.set_yticks(np.arange(len(growth_rates_list))-0.5)
ax4_twin.set_yticklabels(growth_rates_list)
ax4_twin.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
ax4_twin.set_ylim(0,len(df_mem_.condition.unique()))
ax4_twin.xaxis.set_tick_params(labelsize=5)
ax4_twin.yaxis.set_tick_params(labelsize=5)

plt.tight_layout()
fig.savefig('../../figures/fig5abc_SV_scaling_COG_characterization.pdf')


# %%
