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

# Compute the continuum estimate
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
m_aa = 110 / 6E11
atp_frac = 0.8
atp_per_aa = 4
prot_per_atp = 4
r_atp = 300
r_etc = 1500

N_synthase = (cell_mass * theta_dry * theta_prot * atp_per_aa) / (m_aa * atp_frac * r_atp * t_double)
N_respirasome = (cell_mass * theta_dry * theta_prot * atp_per_aa * prot_per_atp) / (m_aa * atp_frac * r_etc * t_double)

# ATP synthesis
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
base_name = 'atp_synthase'
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel(r'F$_1$-F$_0$ ATP synthases per cell')
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.set_ylim([1E2, 1E5])

if ESTIMATE:
    base_name += '_estimate'
    # Plot the continuum
    _ = ax.plot(growth_rate, N_synthase, '-', lw=3, color='grey', alpha=0.4,
                label='scaling with cell size')

    # Plot the point estimate
    _ = ax.hlines(3000, 0, 0.5,'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.vlines(0.5, 0, 3000, 'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.plot(0.5, 3000, 'o', color=colors['dark_brown'], alpha=0.5,
                label='point estimate', ms=6)


if DATA:
    base_name += '_data'
    for g, d in data[data['shorthand']=='atp_synthase'].groupby(['dataset', 'dataset_name']):
        _ = ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=6, markeredgewidth=0.5,
        markeredgecolor='k', color=dataset_colors[g[0]], label=g[1])
    

_ = ax.legend()
plt.savefig(f'../../../figures/presentation/{base_name}.pdf', bbox_inches='tight')


# Electron Transport Complexes
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
base_name = 'respiration'
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel(r'electron transport chains per cell')
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.set_ylim([1E2, 1E5])

if ESTIMATE:
    base_name += '_estimate'
    # Plot the continuum
    _ = ax.plot(growth_rate, N_respirasome, '-', lw=3, color='grey', alpha=0.4,
                label='scaling with cell size')

    # Plot the point estimate
    _ = ax.hlines(2500, 0, 0.5,'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.vlines(0.5, 0, 2500, 'k', linestyle='--', lw=1, label='__nolegend__')
    _ = ax.plot(0.5, 2500, 'o', color=colors['dark_brown'], alpha=0.5,
                label='point estimate', ms=6)


if DATA:
    base_name += '_data'
    for g, d in data[data['shorthand']=='proton_gradient'].groupby(['dataset', 'dataset_name']):
        _ = ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=6, markeredgewidth=0.5,
        markeredgecolor='k', color=dataset_colors[g[0]], label=g[1])
    

_ = ax.legend()
plt.savefig(f'../../../figures/presentation/{base_name}.pdf', bbox_inches='tight')


# %%

# %%
