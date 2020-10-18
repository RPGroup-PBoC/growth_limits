#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size
import scipy.optimize
import scipy.stats
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
merged = tot_prot.merge(ribos)

# Do a simple linear regression on log transform
naa_popt = scipy.stats.linregress(np.log(merged['n_ribosomes'].values),
                              np.log(merged['naa'].values))
vol_popt = scipy.stats.linregress(np.log(merged['n_ribosomes'].values),
                              np.log(merged['volume'].values))

# Results are
naa_slope = naa_popt.slope # Slope is 0.4905
naa_intercept = naa_popt.intercept # Intercept is 16.107
vol_slope = vol_popt.slope
vol_intercept = vol_popt.intercept

# %%
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
for a in ax:
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim([5E3, 5E5])
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# ax.set_ylim([1E8, 5E9])
for g, d in merged.groupby(['dataset', 'dataset_name']):
    ax[0].plot(d['n_ribosomes'],  d['naa'], 'o', color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1], ms=4.5)
    ax[1].plot(d['n_ribosomes'],  d['volume'], 'o', color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1], ms=4.5)
ribosomes = np.logspace(2, 6,500)
naa = ribosomes**naa_slope * np.exp(naa_intercept)
vol = ribosomes**vol_slope * np.exp(vol_intercept)
ax[0].plot(ribosomes, naa, 'k--', lw=1, label='$R^{%s}e^{%s}$' %(np.round(naa_slope, decimals=1),
                                                                  np.round(naa_intercept, decimals=1)))
ax[1].plot(ribosomes, vol, 'k--', lw=1, label='$R^{%s}e^{%s}$' %(np.round(vol_slope, decimals=1),
                                                                  np.round(vol_intercept, decimals=1)))
ax[0].set_xlabel('ribosomes per cell', fontsize=6)
ax[0].set_ylabel('number of amino acids', fontsize=6)
ax[1].set_xlabel('ribosomes per cell', fontsize=6)
ax[1].set_ylabel('cell volume [µm$^3$]', fontsize=6)
ax[0].legend(fontsize=6, loc='lower right')
ax[1].legend(fontsize=6, loc='lower right')
plt.tight_layout()
plt.savefig('../../figures/figA5_linear_regression_naa_volume_ribosomes.pdf')

# %%
