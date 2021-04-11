#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# Load the data and restrict.
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
lipid = data[data['shorthand']=='fas']
pg = data[data['shorthand']=='transpeptidases']

# Compute the scaling relations.
growth_rate = constants['growth_rate']['value']
surface_area = constants['surface_area']['value']
t_double = constants['t_double']['value']

rho_pg = 0.25 # in pg per fL
A_lipid = 0.5 / 1E6 # in square microns
N_leaflet = 4
kcat_acp = 1
kcat_tpd = 2
xlink_frac = 0.2
w_pg = 0.005 # thickness of pg in um
m_murein = 1000 / 6E11 # Mass of murein monomer in pg
theta_lipid = 0.4
N_fabs = (theta_lipid * surface_area * N_leaflet) / (A_lipid * kcat_acp * t_double)
N_tpds = (xlink_frac * w_pg * surface_area * rho_pg) / (m_murein * kcat_tpd * t_double)

# Generate the figures
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.set_xlim([0, 2])
ax.set_yscale('log')
ax.set_ylim([1E2, 1E5])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ACP dehydratases per cell\n(FabZ + FabA)', fontsize=6)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)

# Plot the scaling relationship
ax.plot(growth_rate, N_fabs, '-', lw=3, color='grey', label='surface area scaling',
        alpha=0.4)

# Plot the prediction
ax.plot(0.5, 4E3, 'o', ms=5,  color=colors['dark_brown'],  label='point estimate',
        alpha=0.4)
ax.hlines(4E3, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax.vlines(0.5, 0, 4E3, 'k', linestyle='--', lw=0.75, label='__nolegend__')

# plot the data
for g, d in lipid.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4,
            color=dataset_colors[g[0]], alpha=0.75, markeredgecolor='k',
            markeredgewidth=0.5, label=g[1])
ax.legend(fontsize=6, loc='lower right')
plt.savefig('../../figures/fig3a_lipid_synthesis_plots.svg', bbox_inches='tight')
# %%
# Generate the figures
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.set_xlim([0, 2])
ax.set_yscale('log')
ax.set_ylim([5E0, 5E3])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('murein transpeptidases per cell', fontsize=6)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)

# Plot the scaling relationship
ax.plot(growth_rate, N_tpds, '-', lw=3, color='grey', label='surface area scaling',
        alpha=0.4)

# Plot the prediction
ax.plot(0.5, 100, 'o', ms=5,  color=colors['dark_brown'],  label='point estimate',
        alpha=0.4)
ax.hlines(100, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax.vlines(0.5, 0, 100, 'k', linestyle='--', lw=0.75, label='__nolegend__')

# plot the data
for g, d in pg.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4,
            color=dataset_colors[g[0]], alpha=0.75, markeredgecolor='k',
            markeredgewidth=0.5, label=g[1])
ax.legend(fontsize=6, loc='upper left', ncol=2)
plt.savefig('../../figures/fig3b_pg_biosynthesis.svg', bbox_inches='tight')



# %%
