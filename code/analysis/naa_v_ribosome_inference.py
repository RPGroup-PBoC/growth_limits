#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import prot.viz
import prot.size
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()

# Load data for ribosomes
_data = pd.read_csv('../../data/compiled_estimate_categories.csv')
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
tot_prot = tot_prot[['dataset', 'dataset_name', 'growth_rate_hr', 'condition', 'naa']]

# Combine the datasets
merged = tot_prot.merge(ribos)

# %%
fig, ax = plt.subplots(1, 1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([5E3, 1E5])
ax.set_ylim([1E7, 1E10])
for g, d in merged.groupby(['dataset', 'dataset_name']):
    ax.plot(d['n_ribosomes'], d['naa'], 'o', color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])
# %%
