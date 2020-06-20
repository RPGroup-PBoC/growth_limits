# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size
import prot.estimate
constants = prot.estimate.load_constants()
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')
data = data[data['shorthand']=='sulfur_tport']

# Define constants
theta_S = constants['dry_mass_frac']['value'] * constants['theta_S']['value']
rho = constants['density']['value']
vol = constants['volume']['value']
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
mass = constants['cell_mass']['value']
m_sulf = 32/6E11 # in pg
r_sulf = 10 # in C / s
N_tporters = (theta_S * mass)/ (m_sulf * r_sulf * t_double)

# Set up the figure canvas. 
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([1E1, 5E4])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('number of CysUWA sulfate\ntransporters per cell', fontsize=6)

# Plot the scaling argument
ax.plot(0.5, 1E3, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax.vlines(0.5, 1E1, 1E3, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax.hlines(1E3, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.plot(growth_rate, N_tporters, '-', lw=3, color='grey', label='cell size dependence',
alpha=0.3)
ax.legend(ncol=2, fontsize=6)
plt.savefig('../../figures/fig3b_sulf_transporters.svg', bbox_inches='tight')
# %%
