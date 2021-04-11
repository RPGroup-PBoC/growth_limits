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
fig = plt.figure()
widths = [5, 6]
# heights = [2, 0.5, 0.5, 0.25]
# widths = [5, 5]
heights = [1.75, 1.25, 1.5, 1.5]
spec = fig.add_gridspec(ncols=2, nrows=4, width_ratios=widths)


ax1 = fig.add_subplot(spec[:2,0])
ax2 = fig.add_subplot(spec[2:,0])


# Load the data set
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# total fg per cell cytosol GO:0005829 - cytosol
data_cytosol = data[data.go_terms.str.contains('GO:0005829')]

data_cytosol = data_cytosol.replace('Not Assigned', 'poorly characterized or not assigned')
data_cytosol = data_cytosol.replace('poorly characterized', 'poorly characterized or not assigned')

cog_class_order = ['metabolism',
                'information storage and processing',
                'cellular processes and signaling',
                'poorly characterized or not assigned']
order_dict = dict(zip(cog_class_order,
                np.arange(4)))

color_dict = dict(zip(cog_class_order,
                    ['#679B48', '#D3B15E', '#BF703A', '#788FBD']))

marker_dict = dict(zip(data.dataset.unique(),
                    ['o', 's', 'd', 'v']))

label_dict = dict(zip(data.dataset_name.unique(),
                    ['o', 's', 'd', 'v']))
######################
# relative dist of different COG categories: Plot (D)
######################


# Compute the mass fraction
mass_frac = []
for g, d in data_cytosol.groupby(['dataset', 'condition', 'growth_rate_hr']):
    tot_mass = d['fg_per_cell'].sum()
    sector_mass = d.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr', 'cog_class'])['fg_per_cell'].sum().reset_index()
    frac = sector_mass['fg_per_cell'].values / tot_mass
    sector_mass['frac'] = frac
    sector_mass['dataset_name'] = sector_mass['dataset_name']
    mass_frac.append(sector_mass)
mass_frac = pd.concat(mass_frac, sort=False)

for g, d in mass_frac.groupby(['dataset', 'cog_class', 'dataset_name']):
    ax1.plot(d.growth_rate_hr, d['frac'], marker = marker_dict[g[0]], color=color_dict[g[1]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10, linewidth = 0)

ax1.set_ylim(0,1)
ax1.set_xlim(0,2)
ax1.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax1.set_ylabel('cytosolic protein mass\n fraction (GO term : 0005829)', fontsize=6)
ax1.xaxis.set_tick_params(labelsize=5)
ax1.yaxis.set_tick_params(labelsize=5)

legend_elements = [Line2D([0], [0],
                    marker=label_dict[i],
                    color='w', linewidth = 0,
                    label=i, markeredgecolor='k',
                    markeredgewidth=0.25,
                    markerfacecolor='gray',
                    markersize=4) for i in data.dataset_name.unique()]

ax1.legend(handles=legend_elements, loc='upper left', fontsize = 6)


######################
# metabolism vs information processing; Plot (E)
######################

cm = plt.cm.viridis

for g, d in mass_frac.groupby(['dataset', 'growth_rate_hr', 'dataset_name', 'condition']):
    ax2.scatter(d[d.cog_class == 'information storage and processing']['frac'],
                    d[d.cog_class == 'metabolism']['frac'], marker = marker_dict[g[0]],
                    c = [float(i) for i in d.growth_rate_hr.unique()], cmap = cm,
                    # c = float(d.growth_rate_hr.unique()), cmap=cm,#'viridis',
                    vmin=0, vmax=2, alpha=0.75, edgecolors='k', linewidths=0.25)

ax2.set_ylim(0.25,0.65)
ax2.set_xlim(0.25,0.55)
ax2.set_xlabel('proteomic mass fraction,\ninformation storage and processing', fontsize=6)
ax2.set_ylabel('proteomic mass fraction,\nmetabolism', fontsize=6)
ax2.xaxis.set_tick_params(labelsize=5)
ax2.yaxis.set_tick_params(labelsize=5)

# add in colorbar for growth rate, 0-2 hr-1
cbaxes = inset_axes(ax2, width="30%", height="3%", loc=3)
norm = mpl.colors.Normalize(vmin=0,vmax=2)
sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[0.,2], orientation='horizontal')
cbaxes.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
cbaxes.xaxis.set_tick_params(labelsize=5)
cbaxes.yaxis.set_tick_params(labelsize=5)
cbaxes.xaxis.tick_top()
cbaxes.xaxis.set_label_position('top')
cb.outline.set_linewidth(0)


plt.tight_layout()
fig.savefig('../../figures/fig5de_cytosol.pdf')

# %%
