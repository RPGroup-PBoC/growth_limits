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
# Isolate the glucose transporters (ptsG-crr) and compute mean copy number
GLUCOSE_TRANSPORT_CPLX = ['CPLX-165', 'CPLX-157']
gluc = data[data['complex'].isin(GLUCOSE_TRANSPORT_CPLX)]
gluc = gluc.groupby(['dataset_name', 'growth_rate_hr', 
                   'condition', 'color', 'complex'])['n_units'].mean().reset_index()
gluc = gluc.groupby(['dataset_name', 'growth_rate_hr', 
                   'condition', 'color'])['n_units'].sum().reset_index()

# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
for a in [ax]:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_yscale('log')
    a.set_ylim([2E3, 1E5])
    a.set_yticks([3E3, 1E4, 3E4, 1E5])
    a.set_yticklabels([r'$3\times10^3$', r'$1\times 10^4$', r'$3\times10^4$', 
                       r'$1\times 10^5$'])

# Add labels to axis
a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
a.set_ylabel('complexes per cell', fontsize=8)

for g, d in gluc.groupby(['dataset_name', 'color']):
    _ = ax.plot(d['growth_rate_hr'], d['n_units'], 'o', color=g[1],
               markeredgecolor='k', markeredgewidth=0.35,
               alpha=0.6, label=g[0], ms=5)
leg = ax.legend(fontsize=6)

plt.savefig('../../figures/glucose_transporters.svg')
# %%
WATER_TRANSPORT_CPLX = ['aqpZ', 'CPLX0-7654', 'CPLX0-7534']
water = data[(data['complex'].isin(WATER_TRANSPORT_CPLX)) | (data['gene_name'].isin(WATER_TRANSPORT_CPLX))]
water = water.groupby(['dataset_name', 'growth_rate_hr', 
                   'condition', 'color', 'complex'])['n_units'].mean().reset_index()
water = water.groupby(['dataset_name', 'growth_rate_hr', 
                   'condition', 'color'])['n_units'].sum().reset_index()



# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
for a in [ax]:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_yscale('log')
    a.set_ylim([1E2, 1E5])


# Add labels to axis
a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
a.set_ylabel('complexes per cell', fontsize=8)

for g, d in water.groupby(['dataset_name', 'color']):
    _ = ax.plot(d['growth_rate_hr'], d['n_units'], 'o', color=g[1],
               markeredgecolor='k', markeredgewidth=0.35,
               alpha=0.6, label=g[0], ms=5)
leg = ax.legend(fontsize=6)

plt.savefig('../../figures/water_transporters.svg')
# %%
