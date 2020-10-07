#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

# Load the data and restrict
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
rnap = data[data['shorthand']=='rnap']
sig70 = data[data['shorthand']=='sigma70']


# Compute the predicted trends..
growth_rate = constants['growth_rate']['value']
N_ori = constants['N_ori']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
t_double = constants['t_double']['value']

m_aa = 110 / 6E11 # in pg
m_prot = (300 * 110) / 6E11 # in pg
n_prot = (theta_prot * theta_dry * cell_mass) /m_prot
gamma_mRNA = 1 / 300
n_mRNA = (n_prot / 1E3) * gamma_mRNA
L_mRNA = 1000
L_rRNA = 4500
r_txn = 40 # in nt/s
n_operon = 7
footprint = 80 # in nt
n_tRNA = theta_prot * theta_dry * cell_mass / (m_aa * t_double)
L_tRNA = 80

N_polymerase = (L_rRNA * n_operon  * N_ori/ footprint) + (n_mRNA * L_mRNA / r_txn) + (n_tRNA * L_tRNA) / (r_txn * t_double)


#%%
fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.5))
for a in ax:
    a.plot(growth_rate, N_polymerase, '-', lw=3, color='grey', alpha=0.4, label='replication fork scaling')
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
    a.set_yscale('log')
    a.set_xlim([0, 2])
ax[0].set_ylabel('RNA Polymerases per cell', fontsize=8)
ax[1].set_ylabel('$\sigma^{70}$ (RpoD) per cell', fontsize=8)
ax[0].set_yticks([1E2, 1E3, 1E4, 1E5])
ax[0].set_ylim([1E2, 1E5])
ax[1].set_yticks([1E2, 1E3, 1E4])
ax[1].set_ylim([1E2, 1E4])

# Plot the predictions
for a in ax:
    a.plot(0.5, 1E3, 'o', ms=6, alpha=0.4, color=colors['dark_brown'], label='estimated value')
    a.vlines(0.5, 1, 1E3, color='k', linestyle='--', label='__nolegend__', lw=0.75)
    a.hlines(1E3, 0, 0.5, color='k', linestyle='--', label='__nolegend__', lw=0.75)

# plot the data
for p, a in zip([rnap, sig70], ax.ravel()):
    for g, d in p.groupby(['dataset', 'dataset_name']):
        a.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
                markeredgewidth=0.5, markeredgecolor='k', label=g[1])

# Add legends.
for a in ax:
    a.legend(fontsize=6, loc='lower right')
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

plt.tight_layout()
plt.savefig('../../figures/fig7_RNA_synthesis_plots.svg', bbox_inches='tight')
# %%
