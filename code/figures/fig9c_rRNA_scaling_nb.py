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
n_rRNA_full = n_operon * k_rrna * n_ori * t_double
n_rRNA_half = 0.5*n_operon * k_rrna * n_ori * t_double
n_rRNA_single = k_rrna * n_ori * t_double
n_rRNA_limit = n_operon * k_rrna  * t_double

fig, ax = plt.subplots(1,  1, figsize=(5, 4))
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel(r'number of ribosomal units')
ax.set_ylim([5E3, 8E5])
# ax2 = ax.twinx()
ax.plot(growth_rate, n_rRNA_full, '-', lw=1, color=colors['blue'],
        label=r'rRNA production limit (parallelized, 7 operon copies)')

for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o',  ms=4, markeredgecolor='k',
            markeredgewidth=0.5, color=dataset_colors[g[0]], label=g[1])


# ax.legend(fontsize=6, bbox_to_anchor=(1, 1))
plt.savefig('../../figures/fig9c_rRNA_scaling.pdf')

# %%
