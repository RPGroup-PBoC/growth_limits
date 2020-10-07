#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}

_carbon = data[data['shorthand']=='carbon_tport']
_sulfur = data[data['shorthand']=='sulfur_tport']
_phospho = data[data['shorthand']=='phosphate_tport']




# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([1E2, 1E5])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('mean number of carbohydrate\n transporters per cell', fontsize=6)

# Plot the prediction.
ax.plot(0.5, 1.5E3, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax.vlines(0.5, 1E2, 1.5E3, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax.hlines(1.5E3, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in _carbon.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.legend(ncol=2, fontsize=6)
plt.savefig('../../figures/carbohydrate_tporters_plots.svg', bbox_inches='tight')
# %%


fig, ax = plt.subplots(1, 1, figsize=(3, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 2])
# ax.set_ylim([1E2, 1E5])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('sulfate transporters per cell', fontsize=6)

# Plot the prediction.
ax.plot(0.5, 1E3, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax.vlines(0.5, 1E1, 1E3, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax.hlines(1E3, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in _sulfur.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.legend(ncol=2, fontsize=6)
plt.savefig('../../figures/sulfate_tporters_plots.svg', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([5E1, 1E4])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('phosphate transporters per cell', fontsize=6)

# Plot the prediction.
ax.plot(0.5, 150, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax.vlines(0.5, 1E1, 150, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax.hlines(150, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in _phospho.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.legend(ncol=2, fontsize=6)
plt.savefig('../../figures/phosphate_transporters_plot.svg', bbox_inches='tight')
# %%
# %%
