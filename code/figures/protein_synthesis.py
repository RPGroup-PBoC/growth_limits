#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# load the data 
data = pd.read_csv('../../data/compiled_estimate_categories.csv')

dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}



# tRNA SYNTHETASES
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('total number of\ntRNA synthetases per cell', fontsize=6)
ax.set_yscale('log')
ax.set_ylim([5E3, 5E5])
ax.set_xlim([0, 2])
_trna = data[data['shorthand']=='trna']
ax.plot(0.5, 2.5E4,'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='estimated value')
for g, d in _trna.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, alpha=0.75,
            color=dataset_colors[g[0]], markeredgecolor='k', markeredgewidth=0.5,
            label=g[1])

# Plot the dashed lines.
ax.vlines(0.5, 1, 2.5E4, 'k', linestyle='--', label='__nolegend__', lw=0.75)
ax.hlines(2.5E4, 0, 0.5, 'k', linestyle='--', label='__nolegend__', lw=0.75)
ax.legend(fontsize=6)
plt.savefig('../../figures/tRNA_synthase_plot.svg', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribosomes per cell', fontsize=6)
ax.set_yscale('log')
ax.set_ylim([1E3, 5E5])
ax.set_xlim([0, 2])
_ribo = data[data['shorthand']=='ribosome']
ax.plot(0.5, 1.5E4, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='estimated value')
for g, d in _ribo.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, alpha=0.75,
            color=dataset_colors[g[0]], markeredgecolor='k', markeredgewidth=0.5,
            label=g[1])

# Plot the dashed lines.
ax.vlines(0.5, 1, 1.5E4, 'k', linestyle='--', label='__nolegend__', lw=0.75)
ax.hlines(1.5E4, 0, 0.5, 'k', linestyle='--', label='__nolegend__', lw=0.75)
ax.legend(fontsize=6)
plt.savefig('../../figures/ribosome_plot.svg', bbox_inches='tight')

# %%
