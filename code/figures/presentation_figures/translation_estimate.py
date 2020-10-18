#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

# Define the plotting constants. 
ESTIMATE = True
DATA = False
data = pd.read_csv('../../../data/compiled_estimate_categories.csv')


# Compute the continuum estimates. 
growth_rate = constants['growth_rate']['value'] 
t_double = constants['t_double']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
m_aa = 110 / 6E11
r_translation = 15 # in AA per sec
r_charging = 20
N_ribosomes = (cell_mass * theta_dry * theta_prot) / (m_aa * r_translation * t_double)
N_synthase = (cell_mass * theta_dry * theta_prot) / (m_aa * r_charging * t_double)


fig, ax = plt.subplots(1, 1, figsize=(6, 4))
base_name = 'tRNA_charging'
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('number of tRNA synthetases')
ax.set_xlim([0, 2])
ax.set_ylim([1E3, 5E5])
ax.set_yscale('log')

if ESTIMATE:
    base_name += '_estimate'
    _ = ax.plot(growth_rate, N_synthase, '-', lw=3, color='grey', alpha=0.4,
                label='scaling with cell size')
    _ = ax.plot(0.5, 1E4, 'o', ms=6, color=colors['dark_brown'], alpha=0.4,
                label='point estimate')
    _ = ax.vlines(0.5, 0, 1E4, 'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.hlines(1E4, 0, 0.5, 'k', linestyle='--', lw=1, label='__nolegend__')

if DATA:
    base_name += '_data'
    for g, d in data[data['shorthand']=='trna'].groupby(['dataset', 'dataset_name']):
        _ = ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=6,
                    color=dataset_colors[g[0]], markeredgewidth=0.5,
                    markeredgecolor='k', label=g[1])

_ = ax.legend()
plt.savefig(f'../../../figures/presentation/{base_name}.pdf', bbox_inches='tight')

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
base_name = 'ribosomes'
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('number of ribosomes')
ax.set_xlim([0, 2])
ax.set_ylim([1E3, 5E5])
ax.set_yscale('log')

if ESTIMATE:
    base_name += '_estimate'
    _ = ax.plot(growth_rate, N_synthase, '-', lw=3, color='grey', alpha=0.4,
                label='scaling with cell size')
    _ = ax.plot(0.5, 1E4, 'o', ms=6, color=colors['dark_brown'], alpha=0.4,
                label='point estimate')
    _ = ax.vlines(0.5, 0, 1E4, 'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.hlines(1E4, 0, 0.5, 'k', linestyle='--', lw=1, label='__nolegend__')

if DATA:
    base_name += '_data'
    for g, d in data[data['shorthand']=='ribosome'].groupby(['dataset', 'dataset_name']):
        _ = ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=6,
                    color=dataset_colors[g[0]], markeredgewidth=0.5,
                    markeredgecolor='k', label=g[1])

_ = ax.legend()

plt.savefig(f'../../../figures/presentation/{base_name}.pdf', bbox_inches='tight')

# %%
