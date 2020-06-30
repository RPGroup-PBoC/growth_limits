#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()

# Load the data and restrict
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
dnap = data[data['shorthand']=='dnap']
rnr = data[data['shorthand']=='dntp']

# Compute the cell size dependence. 
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
n_ori = constants['N_ori']['value']
L_genome = 4.6E6 # in nt
r_rnr = 10
r_dna = 600 # in nt per sec
n_pol = 2 # per replication fork
n_fork = 2
N_rnr = 2 * n_ori * L_genome  / (r_rnr * t_double)
N_dnap = n_fork * n_pol * n_ori  
# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribonucleotide reductases per cell', fontsize=6)

# Format the axes
ax.set_yscale('log')
ax.set_ylim([1E1, 1E4])
ax.set_xlim([0, 2])

ax.plot(growth_rate, N_rnr,'-', color='grey', lw=3, alpha=0.5, label='replication fork dependence')

# Plot the predictions
ax.plot(0.5, 200, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax.hlines(200, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax.vlines(0.5, 10, 200, 'k', linestyle='--', lw=0.75, label='__nolegend__')

for g, d in rnr.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.legend(fontsize=6, ncol=2)
plt.savefig('../../figures/fig6a_dNTP_plots.svg', bbox_inches='tight')

# %%



# %%
fig, ax = plt.subplots(2, 1, figsize=(3, 3))
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('DNA polymerase III\nper cell', fontsize=6)
ax[1].set_ylabel('DNA polymerase III\nconcentration [nM]', fontsize=6)

# Format the axes
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_ylim([1, 5E2])
ax[0].set_xlim([0, 2])
ax[1].set_ylim([1, 500])
ax[1].set_xlim([0, 2])


# Plot the predictions
ax[0].plot(growth_rate, N_dnap, '-', lw=3, color='grey', alpha=0.3, label='replication fork dependence')
ax[0].plot(0.5, 6.5, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax[0].hlines(6, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax[0].vlines(0.5, 1, 6.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')



# Plot the concentration range of 50 - 200 nM
ax[1].fill_between([0, 2], 50, 200, color='k', alpha=0.25, label='DNA Pol. III H.E. binding affinity\n (Ason et al. 2000)')

for g, d in dnap.groupby(['dataset', 'dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])
    ax[1].plot(d['growth_rate_hr'], d['concentration_uM'] * 1E3, 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label='__nolegend__')
ax[1].legend(fontsize=6)
ax[0].legend(fontsize=6, ncol=2)
plt.tight_layout()
plt.savefig('../../figures/fig6bc_DNA_polymerase_plots.svg', bbox_inches='tight')
# 



# %%
