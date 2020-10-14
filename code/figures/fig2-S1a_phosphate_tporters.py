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
data = data[data['shorthand']=='phosphate_tport']

# Define constants
theta_P = constants['dry_mass_frac']['value'] * constants['theta_P']['value']
rho = constants['density']['value']
vol = constants['volume']['value']
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
mass = constants['cell_mass']['value']
m_phos = 30/5E11 # in pg
r_phos = 300 # in C / s
N_tporters = (theta_P * mass)/ (m_phos * r_phos * t_double)


# Set up the figure canvas. 
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([1E1, 1E4])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('number of PitA + PitB phosphate\ntransporters per cell', fontsize=6)

# Plot the scaling argument
ax.plot(0.5, 2E2, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='point estimate',
        zorder=1000)
ax.vlines(0.5, 1E1, 2E2, color='k', linestyle='--', lw=0.75, label='__nolegend__',
        zorder=999)
ax.hlines(2E2, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__', 
        zorder=999)

# Plot the data
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.plot(growth_rate, N_tporters, '-', lw=3, color='grey', label='cell size dependence',
alpha=0.3)
ax.legend(ncol=2, fontsize=6)
plt.savefig('../../figures/fig2-S1a_phos_transporters.svg', bbox_inches='tight')

# %%
