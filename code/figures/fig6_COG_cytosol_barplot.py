#%%
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import prot.viz
import prot.size as size
colors, palette = prot.viz.bokeh_theme()
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
prot.viz.plotting_style()


######################
# plot configuration #
######################
fig = plt.figure()#figsize = (6,5))#constrained_layout=True)
widths = [5, 6]
# heights = [2, 0.5, 0.5, 0.25]
# widths = [5, 5]
heights = [1.75, 1.25, 1.5, 1.5]
spec = fig.add_gridspec(ncols=2, nrows=4, width_ratios=widths)
                          # height_ratios=heights)


ax1 = fig.add_subplot(spec[:2,0])
ax2 = fig.add_subplot(spec[2:,0])

ax3 = fig.add_subplot(spec[0,1])
ax4 = fig.add_subplot(spec[1:,1])


# Load the data set
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# total fg per cell cytosol GO:0005829 - cytosol
data_cytosol = data[data.go_terms.str.contains('GO:0005829')]

data_cytosol = data_cytosol.replace('Not Assigned', 'poorly characterized or not assigned')
data_cytosol = data_cytosol.replace('poorly characterized', 'poorly characterized or not assigned')
#

# cog_class_order = ['metabolism',
#                 'cellular processes and signaling',
#                 'information storage and processing',
#                 'poorly characterized or not assigned']
cog_class_order = ['metabolism',
                'information storage and processing',
                'cellular processes and signaling',
                'poorly characterized or not assigned']
order_dict = dict(zip(cog_class_order,
                np.arange(4)))

color_dict = dict(zip(cog_class_order,
                    ['#679B48', '#D3B15E', '#BF703A', '#788FBD']))
                    # ['#679B48', '#BF703A', '#D3B15E', '#788FBD']))


marker_dict = dict(zip(data.dataset.unique(),
                    ['o', 's', 'd', 'v']))

label_dict = dict(zip(data.dataset_name.unique(),
                    ['o', 's', 'd', 'v']))
# ######################
# # Plot 1
# ######################
# # relative dist of different COG categories
#
#
# # Compute the mass fraction
# mass_frac = []
# for g, d in data_cytosol.groupby(['dataset', 'condition', 'growth_rate_hr']):
#     tot_mass = d['fg_per_cell'].sum()
#     sector_mass = d.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr', 'cog_class'])['fg_per_cell'].sum().reset_index()
#     frac = sector_mass['fg_per_cell'].values / tot_mass
#     sector_mass['frac'] = frac
#     sector_mass['dataset_name'] = sector_mass['dataset_name']
#     mass_frac.append(sector_mass)
# mass_frac = pd.concat(mass_frac, sort=False)
#
# for g, d in mass_frac.groupby(['dataset', 'cog_class', 'dataset_name']):
#     ax1.plot(d.growth_rate_hr, d['frac'], marker = marker_dict[g[0]], color=color_dict[g[1]],
#                     alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
#                     label = g[2], ms=4, zorder=10, linewidth = 0)
#
# ax1.set_ylim(0,1)
# ax1.set_xlim(0,2)
# ax1.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
# ax1.set_ylabel('relative cytosolic protein\nfraction (GO term : 0005829)', fontsize=6)
# ax1.xaxis.set_tick_params(labelsize=5)
# ax1.yaxis.set_tick_params(labelsize=5)
#
# legend_elements = [Line2D([0], [0],
#                     marker=label_dict[i],
#                     color='w', linewidth = 0,
#                     label=i, markeredgecolor='k',
#                     markeredgewidth=0.25,
#                     markerfacecolor='gray',
#                     markersize=4) for i in data.dataset_name.unique()]
#
#
#
# ax1.legend(handles=legend_elements, loc='upper left', fontsize = 6)
#
#
# ######################
# # Plot 2
# ######################
# # metabolism vs information processing
#
# cm = plt.cm.get_cmap('viridis')
#
# for g, d in mass_frac.groupby(['dataset', 'growth_rate_hr', 'dataset_name', 'condition']):
#     ax2.scatter(d[d.cog_class == 'information storage and processing']['frac'],
#                     d[d.cog_class == 'metabolism']['frac'], marker = marker_dict[g[0]],
#                     c = d.growth_rate_hr.unique(), cmap=cm,#'viridis',
#                     vmin=0, vmax=2, alpha=0.75, edgecolors='k', linewidths=0.25)
#
# ax2.set_ylim(0.25,0.65)
# ax2.set_xlim(0.25,0.55)
# ax2.set_xlabel('relative mass fraction,\ninformation storage and processing', fontsize=6)
# ax2.set_ylabel('relative mass fraction,\nmetabolism', fontsize=6)
# ax2.xaxis.set_tick_params(labelsize=5)
# ax2.yaxis.set_tick_params(labelsize=5)
#
# cbaxes = inset_axes(ax2, width="30%", height="5%", loc=3)
# cb = plt.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=2), cmap=cm),
#         cax=cbaxes, ticks=[0.,2], orientation='horizontal')
# # cbaxes.ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)#, rotation=270)
# cbaxes.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
# cbaxes.xaxis.set_tick_params(labelsize=5)
# cbaxes.yaxis.set_tick_params(labelsize=5)
# cbaxes.xaxis.tick_top()
# cbaxes.xaxis.set_label_position('top')
# cb.outline.set_linewidth(0)
#
# ######################
# # Plot 3 - protein concentration in cell
# ######################
#
# data_cytosol_fg_summary = pd.DataFrame()
# for c, d in data_cytosol.groupby(['dataset', 'condition', 'growth_rate_hr']):
#     d_ = d[d.gene_name != 'tufA']
#     d_ = d_[d_.gene_name != 'tufB']
#
#     # V = 0.28 * np.exp(1.33 * c[2])
#     # SA = 2 * np.pi *  V**(2/3)
#     w = size.lambda2width(c[2])
#     l = size.lambda2length(c[2])
#     V = size.lambda2size(c[2])
#     SA = size.rod_SA(l, w, V)
#
#     fg_tot = d_.fg_per_cell.sum()
#     fg_SA = fg_tot / V
#
#     data_list = {'fg per um2' : fg_SA,
#                 'growth_rate_hr' : c[2],
#                 'dataset' : c[0],
#                 'condition' : c[1]}
#     data_cytosol_fg_summary = \
#             data_cytosol_fg_summary.append(data_list,
#             ignore_index = True)
#
# for g, d in data_cytosol_fg_summary.groupby(['dataset', 'condition', 'growth_rate_hr']):
#     ax3.plot(g[2], d['fg per um2'], 'o', color=dataset_colors[g[0]],
#                     alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
#                     label = g[0], ms=4, zorder=10)
# ax3.set_ylim(0,210)
# ax3.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
# ax3.set_ylabel('[fg per $\mu m^3$]', fontsize=6)
# ax3.xaxis.set_tick_params(labelsize=5)
# ax3.yaxis.set_tick_params(labelsize=5)
#


######################
# Plot 3
######################
# Plot distribution of inner membrane proteins, GO:0005886
# Load the complex subunit counts.
subunits = pd.read_csv('../../data/compiled_annotated_complexes.csv')

df_cyt_schmidt= data_cytosol#[data_cytosol.dataset == 'schmidt_2016']

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

df_cyt_schmidt_['rel_fg_per_cell'] = df_cyt_schmidt_.groupby(['dataset', 'growth_rate_hr', 'condition']).transform(lambda x: (x / x.sum()))['fg_per_cell']
df_cyt_schmidt_ = df_cyt_schmidt_.sort_values(by=['growth_rate_hr', 'gene_name'], ascending = False)

# y_order = dict(zip(df_cyt_schmidt_.condition.unique(), np.arange(len(df_cyt_schmidt_.condition.unique()))))

gr_yaxis = []
count =  -1
for c, d in df_cyt_schmidt_.groupby(['dataset', 'growth_rate_hr', 'condition'], sort=False):
    gr_yaxis = np.append(gr_yaxis,c[1])
    count += 1
    for c_ in cog_class_order:
        if c_ == 'metabolism':
            resp = d[d.gene_name == 'respiration']
            ax4.barh(count, resp.rel_fg_per_cell, height=0.9, color='#679B48', alpha=0.3, #84A779
                        linewidth=0.1)

            c_uptake = d[d.gene_name == 'carbon_uptake']
            lefts = resp.rel_fg_per_cell.sum()
            ax4.barh(count, c_uptake.rel_fg_per_cell.sum(),  height=0.9, color='#679B48', alpha=0.7,
                            left=lefts, linewidth=0.1)

            # lefts = d[d.gene_name == 'respiration']
            lefts += d[d.gene_name == 'carbon_uptake'].rel_fg_per_cell.sum()
            meta = d[d.gene_name != 'respiration'].copy()
            meta = meta[meta.gene_name !='carbon_uptake'].copy()

            ax4.barh(count, meta.rel_fg_per_cell.sum(), height=0.9, color='#679B48',
                            left=lefts,  linewidth=0.1)

            lefts += meta[meta.cog_class == c_].rel_fg_per_cell.sum()

        else:
            ax4.barh(count, d[d.cog_class == c_].rel_fg_per_cell.sum(), height=0.9,
                        color=color_dict[c_], left=lefts, linewidth=0.1)
            lefts += d[d.cog_class == c_].rel_fg_per_cell.sum()



ax4.set_xlim(0,1)
# ax4.set_ylim(0,count)
ax4.set_ylim(-0.5,count+0.5)
# ax4.set_yticks(np.arange(len(df_cyt_schmidt_.condition.unique()))-0.5)
# ax4.set_yticklabels(df_cyt_schmidt_.condition.unique())
# ax4.set_xlabel('relative cytosolic protein\nabundance (GO term : 0005829)', fontsize=6)
ax4.xaxis.set_tick_params(labelsize=5)
ax4.yaxis.set_tick_params(labelsize=5)

# growth_rates_list = [d.growth_rate_hr.unique()[0] for c, d in df_cyt_schmidt_.groupby('condition', sort=False)]
ax4_twin = ax4.twinx()
# ax4_twin.set_yticks(np.arange(len(growth_rates_list))-0.5)
# ax4_twin.set_yticklabels(growth_rates_list)
ax4_twin.set_yticks(np.arange(len(gr_yaxis))-0.5)
ax4_twin.set_yticklabels(gr_yaxis)
ax4_twin.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
# ax4_twin.set_ylim(0,len(df_cyt_schmidt_.condition.unique()))
ax4_twin.set_ylim(-0.5,count+0.5)
ax4_twin.xaxis.set_tick_params(labelsize=5)
ax4_twin.yaxis.set_tick_params(labelsize=5)
#
plt.tight_layout()
fig.savefig('../../figures/fig6_COG_cytosol_barplot_all.pdf')

# %%
