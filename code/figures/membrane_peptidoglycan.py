#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()


dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
_nitrogen = data[data['shorthand']=='nitrogen_tport']
# Load the data and restrict.  
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
_lipid = data[data['shorthand']=='fas']

fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.set_xlim([0, 2])
ax.set_yscale('log')
ax.set_ylim([5E2, 1E5])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ACP dehydratases per cell\n(FabZ + FabA)', fontsize=6)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)

# Plot the predictions. 
ax.plot(0.5, 3E3, 'o', ms=6,  color=colors['dark_brown'],  label='estimated value',
        alpha=0.75)
ax.hlines(3E3, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax.vlines(0.5,5E2, 3E3, 'k', linestyle='--', lw=0.75, label='__nolegend__')

# plot the data
for g, d in _lipid.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, 
            color=dataset_colors[g[0]], alpha=0.75, markeredgecolor='k',
            markeredgewidth=0.5, label=g[1])

ax.legend(fontsize=6, loc='lower right')
plt.savefig('../../figures/lipid_synthesis_plots.svg', bbox_inches='tight')
# %%
