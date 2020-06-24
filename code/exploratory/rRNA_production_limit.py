#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

# Load experimental data
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
data = data[data['shorthand']=='ribosome']

# Define constants
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
n_ori = constants['N_ori']['value']

# Compute the maximum number of rRNA produced per doubling time
k_rrna = 1 # functional rRNA units per second at steady state
n_operon = 7 # average number of functional ribosomal operons per chromosome
n_rRNA = n_operon * k_rrna * n_ori * t_double
n_rRNA_limit = n_operon * k_rrna  * t_double

# # Compute the maximum number of 
# r_tsl = 15 # in AA / s
# n_aa = 7500 # AA per ribosome 
# n_footprint = 30
# ribo_density = n_aa / n_footprint
# t_ribo = n_aa / r_tsl
# r_ribo = ribo_density / t_ribo
# n_mRNA = 10
# n_ribo = r_ribo * n_mRNA * n_ori * t_double 

# Compute an estimate of the nubmer of ribosomes needed as a fn of growth rate
m_aa = 110 / 6E11
m_prot = constants['cell_mass']['value'] * constants['theta_prot']['value'] * constants['dry_mass_frac']['value']
n_aa = m_prot / m_aa
r_trsl = 8# in AA per second
n_ribo = n_aa / (r_trsl * t_double)

fig, ax = plt.subplots(1,  1)
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel(r'number of ribosomal units')
# ax2 = ax.twinx()
ax.plot(growth_rate, n_rRNA, '-', lw=1, color=colors['blue'], 
        label=r'rRNA production limit (parallelized)')
ax.plot(growth_rate, n_rRNA_limit, '--', lw=1, color=colors['blue'], 
        label=r'rRNA production limit (non-parallelized)')
ax.plot(growth_rate, n_ribo, '-', lw=1, color=colors['green'],
        label='predicted number of ribosomes')
# ax2.plot(growth_rate, n_ori, 'k--')
# ax2.set_ylabel('average number of origins')

for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o',  ms=4, markeredgecolor='k',
            markeredgewidth=0.5, color=dataset_colors[g[0]], label=g[1])

# ax.plot(growth_rate, n_ribo, '-', lw=1, color=colors['green'], 
        # label=r'ribosomal protein production limit')


ax.legend()
plt.savefig('./rRNA_scaling.pdf')

# %%
