#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

SHOW_POINT = True
SHOW_CONT = False
SHOW_DATA = True 
base_name = 'carbon_estimate'

# Load the compiled data. 
data = pd.read_csv('../../../data/compiled_estimate_categories.csv')


# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(6,4))
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('number of glucose transporters per cell')
ax.set_yscale('log')
ax.set_ylim([1E2, 1E5])
ax.set_xlim([0, 2])


if SHOW_POINT:
    base_name += '_point'
    _ = ax.plot(0.5, 1300, 'o', color=colors['dark_brown'], alpha=0.5,
                label='point estimate')
    _ = ax.vlines(0.5, 0, 1300, 'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.hlines(1300, 0, 0.5, 'k', linestyle='--', lw=1, label='__nolegend__')

if SHOW_CONT:
    base_name += '_continuum'
    # Compute the continuum estimate.
    growth_rate = constants['growth_rate']['value']
    cell_mass = constants['cell_mass']['value']
    theta_dry = constants['dry_mass_frac']['value']
    theta_c = constants['theta_C']['value']
    t_double = constants['t_double']['value']
    M_carbon = 12 / 6E11 # in pg
    N_carbon = cell_mass * theta_dry * theta_c / M_carbon
    k_tport = 6 * 200 # in C per sec.
    N_tporters = N_carbon / (k_tport * t_double)
    _ = ax.plot(growth_rate, N_tporters, '-', lw=3, color='grey', alpha=0.4,
                label='scaling with cell size')

# Plot the data
if SHOW_DATA:
    base_name += '_data'
    for g, d in data[data['shorthand']=='glucose_tport'].groupby(['dataset', 'dataset_name']):
        _ = ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=6, 
               color=dataset_colors[g[0]], markeredgecolor='k', markeredgewidth=0.5,
               label=g[1])

# Add legend info
_ = ax.legend()
plt.savefig(f'../../../figures/presentation/{base_name}.pdf', bbox_inches='tight')
# %%
