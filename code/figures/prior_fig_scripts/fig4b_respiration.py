#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

# Load the data set
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
data = data[data['shorthand']=='proton_gradient']

# Compute the scaling trend.
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
m_aa = 110 / 6E11 # in pg
r_atp = 300 # per second per synthase
atp_aa = 5
prot_atp = 4
r_etc = 1500

# N_synthase = (cell_mass * theta_dry * theta_prot * atp_aa) / (m_aa * r_atp * t_double)
tot_prot = prot.size.lambda2P(growth_rate) / 1E3
N_synthase = (tot_prot * atp_aa) / (m_aa * r_atp * t_double / np.log(2))
N_ETC = N_synthase * prot_atp * r_atp / r_etc


# Instantiate andf ormat the axis
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('electron transport complexes\nper cell', fontsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([1E2, 5E4])
# Plot the scaling relationship
ax.plot(growth_rate, N_ETC, '-', lw=3, color='grey', alpha=0.4, label='cell size dependence')

# Plot the estimate value
ax.plot(0.5, 2500, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.5, label='estimated value')
ax.vlines(0.5, 0, 2500, 'k', lw=1, linestyle='--', label='__nolegend__')
ax.hlines(2500, 0, 0.5, 'k', lw=1, linestyle='--', label='__nolegend__')

for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.legend(fontsize=6)
plt.savefig('../../figures/fig4b_respiration.svg')



# %%
