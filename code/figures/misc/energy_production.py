#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
prot.viz.plotting_style()

# Load the data set
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
data = data[(data['shorthand']=='atp_synthase') |
            (data['shorthand']=='proton_gradient')]

# %%
fig, ax = plt.subplots(1, 2, figsize=(6.5, 2 ))
axes = {'atp_synthase':ax[0], 'proton_gradient':ax[1]}

# Populate with data
for g, d in data.groupby(['shorthand', 'dataset', 'dataset_name']):
    axes[g[0]].plot(d['growth_rate_hr'], d['n_complex'], 'o', color=dataset_colors[g[1]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4)

# plot the prediction
ax[0].plot(0.5, 3E3, 'o', ms=6, color=colors['dark_brown'], label='estimated value',
           alpha=0.4)
ax[1].plot(0.5, 2500, 'o', ms=6, color=colors['dark_brown'], label='estimated value',
           alpha=0.4)


# Format the axes
for a in ax:
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
    a.set_yscale('log')
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_ylabel('complexes per cell', fontsize=6)
    a.set_ylim([1E2, 5E4])
    a.set_xlim([0, 2])
    a.legend(fontsize=5, loc='lower right')

# Add locations of estimates
ax[0].hlines(3E3, 0, 0.5, color='k', linestyle='--', lw=0.5)
ax[0].vlines(0.5, -0.1, 3E3, color='k', linestyle='--', lw=0.5)
ax[1].hlines(2500, 0, 0.5, color='k', linestyle='--', lw=0.5)
ax[1].vlines(0.5, -0.1, 2500, color='k', linestyle='--', lw=0.5)

plt.tight_layout()
plt.savefig('../../figures/energy_estimate_plots.svg')
# %%
