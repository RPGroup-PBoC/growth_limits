#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()



dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
# %%
# Load the data sets 
data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# %%
fig, ax = plt.subplots(3, 1, figsize=(4, 5))

for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.set_ylabel('complexes per cell', fontsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

for g, d in data.groupby(['dataset', 'dataset_name']):
    _allcarbo = d[d['shorthand']=='carbon_tport_tot']
    _gluc = d[d['shorthand']=='glucose_tport']
    _xyl = d[d['shorthand']=='xylose_tport']

    ax[0].plot(_allcarbo['growth_rate_hr'], _allcarbo['n_complex'], 'o', ms=4,
            markeredgecolor='k', markeredgewidth=0.5, color=dataset_colors[g[0]],
            label=g[1])
    ax[1].plot(_gluc['growth_rate_hr'], _gluc['n_complex'], 'o', ms=4,
            markeredgecolor='k', markeredgewidth=0.5, color=dataset_colors[g[0]],
            label=g[1])
    ax[2].plot(_xyl['growth_rate_hr'], _xyl['n_complex'], 'o', ms=4,
            markeredgecolor='k', markeredgewidth=0.5, color=dataset_colors[g[0]],
            label=g[1])

prot.viz.titlebox(ax[0], 'carbohydrate transporters (total)', color='k', bgcolor=colors['pale_yellow'], size=6)
prot.viz.titlebox(ax[1], 'glucose and mannose transporters (total)', color='k', bgcolor=colors['pale_yellow'], size=6)
prot.viz.titlebox(ax[2], 'xylose transporters (total)', color='k', bgcolor=colors['pale_yellow'], size=6)

ax[0].set_yscale('log')
ax[0].set_ylim([1E3, 5E5])
ax[1].set_yscale('log')
ax[1].set_ylim([1E3, 1E5])
ax[2].set_yscale('log')
ax[2].set_ylim([1, 1E4])
plt.tight_layout()
plt.savefig('../../figures/carbohydrate_transport_detailed.pdf', bbox_inches='tight')
# %%

