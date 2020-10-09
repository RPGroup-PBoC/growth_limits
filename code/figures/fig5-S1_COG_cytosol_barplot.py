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
fig = plt.figure(figsize = (4,4.5))#constrained_layout=True)
# widths = [6, 7]
# heights = [2, 0.5, 0.5, 0.25]
widths = [3]
heights = [1, 6]
spec = fig.add_gridspec(ncols=1, nrows=2, width_ratios=widths,
                          height_ratios=heights)


ax3 = fig.add_subplot(spec[0])
ax4 = fig.add_subplot(spec[1])

######################
# Plot 2
######################
# Load the data set
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# total fg per cell cytosol GO:0005829 - cytosol
data_cytosol = data[data.go_terms.str.contains('GO:0005737')]#0005829')]

data_cytosol_fg_summary = pd.DataFrame()
for c, d in data_cytosol.groupby(['dataset', 'condition', 'growth_rate_hr']):
    d_ = d[d.gene_name != 'tufA']
    d_ = d_[d_.gene_name != 'tufB']

    # V = 0.28 * np.exp(1.33 * c[2])
    # SA = 2 * np.pi *  V**(2/3)
    w = size.lambda2width(c[2])
    l = size.lambda2length(c[2])
    V = size.lambda2size(c[2])
    SA = size.rod_SA(l, w, V)

    fg_tot = d_.fg_per_cell.sum()
    fg_SA = fg_tot / V

    data_list = {'fg per um2' : fg_SA,
                'growth_rate_hr' : c[2],
                'dataset' : c[0],
                'condition' : c[1]}
    data_cytosol_fg_summary = \
            data_cytosol_fg_summary.append(data_list,
            ignore_index = True)

for g, d in data_cytosol_fg_summary.groupby(['dataset', 'condition', 'growth_rate_hr']):
    ax3.plot(g[2], d['fg per um2'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[0], ms=4, zorder=10)
ax3.set_ylim(0,210)
ax3.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax3.set_ylabel('[fg per $\mu m^3$]', fontsize=6)
ax3.xaxis.set_tick_params(labelsize=5)
ax3.yaxis.set_tick_params(labelsize=5)



######################
# Plot 3
######################
# Plot distribution of inner membrane proteins, GO:0005886
# Load the complex subunit counts.
subunits = pd.read_csv('../../data/compiled_annotated_complexes.csv')

df_cyt_schmidt= data_cytosol[data_cytosol.dataset == 'schmidt_2016']

df_cyt_schmidt = df_cyt_schmidt.replace('Not Assigned', 'poorly characterized or not assigned')
df_cyt_schmidt = df_cyt_schmidt.replace('poorly characterized', 'poorly characterized or not assigned')

df_cyt_schmidt_ = pd.DataFrame()
for c,d in df_cyt_schmidt.groupby(['dataset', 'condition',
                                   'growth_rate_hr', 'gene_name']):
    data_list = {'gene_name' : d.gene_name.unique()[0] ,
                 'cog_class' : d.cog_class.unique()[0] ,
                 'condition' : d.condition.unique()[0] ,
                 'dataset' : d.dataset.unique()[0] ,
                 'dataset_name' : d.dataset_name.unique()[0] ,
                 'fg_per_cell' : d.fg_per_cell.sum(),
                'growth_rate_hr' :  d.growth_rate_hr.unique()[0]}
    df_cyt_schmidt_ = df_cyt_schmidt_.append(data_list,
                                            ignore_index=True)

df_cyt_schmidt_['rel_fg_per_cell'] = df_cyt_schmidt_.groupby('condition').transform(lambda x: (x / x.sum()))['fg_per_cell']
df_cyt_schmidt_ = df_cyt_schmidt_.sort_values(by=['growth_rate_hr', 'gene_name'], ascending = False)

y_order = dict(zip(df_cyt_schmidt_.condition.unique(), np.arange(len(df_cyt_schmidt_.condition.unique()))))
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


for c, d in df_cyt_schmidt_.groupby('condition', sort=False):
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
ax4.set_ylim(0,len(df_cyt_schmidt_.condition.unique()))
ax4.set_yticks(np.arange(len(df_cyt_schmidt_.condition.unique()))-0.5)
ax4.set_yticklabels(df_cyt_schmidt_.condition.unique())
ax4.set_xlabel('relative protein abundance', fontsize=6)
ax4.set_xlabel('relative cytoplasmic protein\nabundance (GO term : 0005829)', fontsize=6)
ax4.xaxis.set_tick_params(labelsize=5)
ax4.yaxis.set_tick_params(labelsize=5)

growth_rates_list = [d.growth_rate_hr.unique()[0] for c, d in df_cyt_schmidt_.groupby('condition', sort=False)]
ax4_twin = ax4.twinx()
ax4_twin.set_yticks(np.arange(len(growth_rates_list))-0.5)
ax4_twin.set_yticklabels(growth_rates_list)
ax4_twin.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
ax4_twin.set_ylim(0,len(df_cyt_schmidt_.condition.unique()))
ax4_twin.xaxis.set_tick_params(labelsize=5)
ax4_twin.yaxis.set_tick_params(labelsize=5)
#
plt.tight_layout()
fig.savefig('../../figures/fig5-S1_COG_cytosol_barplot_2.svg')


# %%
