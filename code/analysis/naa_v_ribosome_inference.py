#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import prot.viz
import prot.size
import scipy.optimize
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()

# Load data for ribosomes
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
ribos = data[data['shorthand']=='ribosome']
ribos['conc'] = ribos['n_complex'].values 
ribos = ribos[['dataset', 'dataset_name', 'growth_rate_hr', 'condition', 'conc']]
ribos.rename(columns={'conc':'n_ribosomes'}, inplace=True)

# Load compiled data for total amino acid content. 
_data = pd.read_csv('../../data/compiled_absolute_measurements.csv')
tot_prot = _data.groupby(['dataset', 'dataset_name', 
                          'growth_rate_hr', 'condition'])['fg_per_cell'].sum().reset_index()
tot_prot['volume'] = prot.size.lambda2size(tot_prot['growth_rate_hr'].values)
tot_prot['naa'] = (tot_prot['fg_per_cell'].values * 1E-3 / (110 / 6E11)) 
tot_prot = tot_prot[['dataset', 'dataset_name', 'growth_rate_hr', 'condition', 'naa', 'volume']]

# Combine the datasets




# %%
fig, ax = plt.subplots(1, 1)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim([5E3, 5E4])
# ax.set_ylim([1E8, 5E9])
for g, d in merged.groupby(['dataset', 'dataset_name']):
    ax.plot(d['n_ribosomes'] / d['volume'],  d['volume'], 'o', color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])
ax.set_xlabel('ribosomes per cubic micron')
# ax.set_ylabel('amino acids per cubic micron')
# %%
