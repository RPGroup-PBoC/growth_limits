#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz 
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

#%%
# Load the dataset. 
data = pd.read_csv('../../data/compiled_annotated_complexes.csv')
counts = pd.read_csv('../../data/compiled_absolute_measurements.csv')

#  Assign the colors as necessary.
author_colors = {'Schmidt et al. 2016': colors['red'],
                 'Peebo et al. 2015': colors['yellow'],
                 'Valgepea et al. 2013': colors['blue'],
                 'Li et al. 2014': colors['dark_green']}
for a, c in author_colors.items():
    data.loc[data['dataset_name']==a, 'color'] = c



# %%
DNAP_CPLX = ['CPLX0-3803']
dnap = data[data['complex'].isin(DNAP_CPLX)]
dnap = dnap.groupby(['dataset_name', 'growth_rate_hr', 
                   'condition', 'color'])['n_units'].sum().reset_index()

# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
for a in [ax]:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_yscale('log')
    a.set_ylim([10, 2.5E3])

# Add labels to axis
a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
a.set_ylabel('complexes per cell', fontsize=8)

for g, d in dnap.groupby(['dataset_name', 'color']):
    _ = ax.plot(d['growth_rate_hr'], d['n_units'], 'o', color=g[1],
               markeredgecolor='k', markeredgewidth=0.35,
               alpha=0.6, label=g[0], ms=5)
leg = ax.legend(fontsize=6)
plt.savefig('../../figures/dna_polymerase_III.svg')

# %%
# Ribonucleotide reductase
RNR_CPLX = ['MONOMER0-2863']
rnr = data[data['complex'].isin(RNR_CPLX)]
rnr = rnr.groupby(['dataset_name', 'growth_rate_hr', 
                   'condition', 'color'])['n_units'].sum().reset_index()

# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
for a in [ax]:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    # a.set_yscale('log')
    # a.set_ylim([10, 2.5E3])

# Add labels to axis
a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
a.set_ylabel('complexes per cell', fontsize=8)

for g, d in rnr.groupby(['dataset_name', 'color']):
    _ = ax.plot(d['growth_rate_hr'], d['n_units'], 'o', color=g[1],
               markeredgecolor='k', markeredgewidth=0.35,
               alpha=0.6, label=g[0], ms=5)
leg = ax.legend(fontsize=6)
plt.savefig('../../figures/ribonucleo_reductase.svg')



# %%
# %%
