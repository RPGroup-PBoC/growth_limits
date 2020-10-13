#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

# Load the data
data = pd.read_csv('../../data/compiled_estimate_categories.csv', comment='#')
data = data[data['shorthand']=='trna']

# Compute the trend
growth_rate = constants['growth_rate']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
t_double = constants['t_double']['value']
m_aa = 110 / 6E11 # in pg
k_trna = 20 # per sec

N_synthase = cell_mass * theta_dry * theta_prot / (m_aa * k_trna * t_double)


# Instantiate the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.set_ylim([1E3, 3E5])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('number of tRNA synthetases', fontsize=6)

# Plot the cell size dependence
ax.plot(growth_rate, N_synthase, lw=3, color='grey', alpha=0.4, label='cell size dependence')

# Plot the point estimate.
estimate = 1E4
ax.plot(0.5, estimate, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax.vlines(0.5, 0, estimate, 'k', linestyle='--', lw=1, label='__nolegend__')
ax.hlines(estimate, 0, 0.5, 'k', linestyle='--', lw=1, label='__nolegend__')

# Plot the experimetnal data
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.legend(fontsize=6)
plt.savefig('../../figures/fig9-S1_tRNA_synthases.svg', bbox_inches='tight')
# %%


# %%
