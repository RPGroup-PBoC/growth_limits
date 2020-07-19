
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')
_carbon = data[data['shorthand']=='carbon_tport']

# Define constants
theta_C = constants['dry_mass_frac']['value'] * constants['theta_C']['value']
rho = constants['density']['value']
vol = constants['volume']['value']
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
mass = constants['cell_mass']['value']
m_carbon = 12/6E11 # in pg
r_carbon = 1000 # in C / s
N_tporters = (theta_C * mass)/ (m_carbon * r_carbon * t_double)

# Set up the figure canvas. 
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([1E2, 5E5])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('total number of \n PTS transporters per cell', fontsize=6)

# Plot the scaling argument
ax.plot(0.5, 2E3, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax.vlines(0.5, 1E2, 2E3, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax.hlines(2E3, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in _carbon.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.plot(growth_rate, N_tporters, '-', lw=3, color='grey', label='cell size dependence',
alpha=0.3)
ax.legend(ncol=2, fontsize=6)
plt.savefig('../../figures/fig2a_carbon_tporters.svg', bbox_inches='tight')


# %%
