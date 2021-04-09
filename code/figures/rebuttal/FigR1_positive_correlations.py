#%%
# from re import 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from seaborn.axisgrid import jointplot 
import prot.viz 
import tqdm 
import scipy.stats
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
#%%
# Load the individual abundance data 
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv', comment='#')


#%%
# FOr each protein, dataset, and condition, compute the total mass fractions.
dfs = []
for g, d in data.groupby(['condition', 'dataset', 'dataset_name', 'growth_rate_hr']):
    tot_mass = d['fg_per_cell'].sum()
    d['mass_fraction'] = d['fg_per_cell'].values / tot_mass 
    dfs.append(d)
mass_fracs = pd.concat(dfs)

# for each gene, compute the correlation coefficient with growth rate 
pearson = pd.DataFrame([])
for g, d in tqdm.tqdm(mass_fracs.groupby(['gene_name', 'dataset', 'dataset_name']),
                      desc='Computing linear correlations'):
    d.dropna(inplace=True)
    if len(d) >= 2:
        r, p = scipy.stats.pearsonr(d['growth_rate_hr'], d['mass_fraction'])
        pearson = pearson.append({
                        'gene_name': g[0],
                        'dataset': g[1],
                        'dataset_name': g[2],
                        'pearson_r': r,
                        'p_value': p},
                        ignore_index=True) 

#%%
# Sort the genes by their correlation 
pearson.dropna(inplace=True)
pearson.sort_values(by='pearson_r', inplace=True, ascending=False)
signif = pearson[pearson['p_value'] <= 0.05]
pos = signif[signif['pearson_r'] > 0.75]

# %%
# Get top 16 genes that are not ribosomal. 
ribo_prots = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']
genes = pos['gene_name'].unique()
genes = [g for g in genes if g not in ribo_prots]
# Get the top 16 non ribosomal proteins
#%%
ind = 80 
top = genes[ind:ind + 16]

pos_data = mass_fracs[mass_fracs['gene_name'].isin(top)]

fig, ax = plt.subplots(4, 4, figsize=(8, 8), sharex=True)
for i in range(4):
    ax[-1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
for a in ax.ravel():
    a.set_ylabel('mass fraction', fontsize=8)
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

iter = 0
for g, d in pos_data.groupby(['gene_name']):
    _ax = ax.ravel()[iter]
    _ax.set_title(g)
    for _g, _d in d.groupby(['dataset', 'dataset_name']):
        _ax.plot(_d['growth_rate_hr'], _d['mass_fraction'], 'o', ms=3, 
                color=dataset_colors[_g[0]], alpha=0.5, markeredgewidth=1,
                markeredgecolor='k')
    iter += 1
plt.tight_layout()
# %%

# Look at specific proteins 
other_prots = ['iadA', 'ychF', 'wecB', 'carB']
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax = ax.ravel()
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_ylabel('mass fraction', fontsize=8)
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)

_data = mass_fracs[mass_fracs['gene_name'].isin(other_prots)]
iter = 0
for g, d in _data.groupby(['gene_name']):
    _ax = ax[iter]
    for _g, _d in d.groupby(['dataset', 'dataset_name']):
        _ax.plot(_d['growth_rate_hr'], _d['mass_fraction'], 'o', ms=3, 
                color=dataset_colors[_g[0]], alpha=0.5, markeredgewidth=1,
                markeredgecolor='k', label=_g[1])
    _ax.set_title(g, fontsize=8)

    iter += 1
ax[0].legend(fontsize=6)
plt.tight_layout()
plt.savefig('../../../FigRX_positive_correlations.pdf', bbox_inches='tight')
# %%
