#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Load the data and restrict
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
rnap = data[data['shorthand']=='rnap']
sig70 = data[data['shorthand']=='sigma70']

dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}

fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.5))
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
    a.set_yscale('log')
    a.set_xlim([0, 2])
ax[0].set_ylabel('RNA Polymerases per cell', fontsize=8)
ax[1].set_ylabel('$\sigma^{70}$ (RpoD) per cell', fontsize=8)
ax[0].set_yticks([1E2, 1E3, 1E4, 1E5])
ax[0].set_ylim([1E2, 1E5])
ax[1].set_yticks([1E2, 1E3, 1E4])
ax[1].set_ylim([1E2, 1E4])

# Plot the predictions
for a in ax:
    a.plot(0.5, 500, 'o', ms=6, color=colors['dark_brown'], label='estimated value')
    a.vlines(0.5, 1, 500, color='k', linestyle='--', label='__nolegend__', lw=0.75)
    a.hlines(500, 0, 0.5, color='k', linestyle='--', label='__nolegend__', lw=0.75)

# plot the data
for p, a in zip([rnap, sig70], ax.ravel()):
    for g, d in p.groupby(['dataset', 'dataset_name']):
        a.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
                markeredgewidth=0.5, markeredgecolor='k', label=g[1])

# Add legends. 
for a in ax:
    a.legend(fontsize=6)
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

plt.tight_layout()
plt.savefig('../../figures/RNA_synthesis_plots.svg', bbox_inches='tight')
# %%
