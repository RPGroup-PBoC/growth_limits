#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/compiled_estimate_categories.csv')


dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}


# %%
# Isolate the necessary data sets. 
rna = data[(data['shorthand']=='rnap') | (data['shorthand']=='all_sigma')]
proteins = data[(data['shorthand']=='ribosome') | (data['shorthand']=='trna')]
atp = data[(data['shorthand']=='atp_synthase') | (data['shorthand']=='proton_gradient')]

# %%
# plot teh RNAP ratio first. 
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('total RNA polymerase / total σ-factor', fontsize=6)
# ax.set_yscale('log')
ax.set_ylim([1, 10])

for g, d in rna.groupby(['dataset', 'dataset_name']):
    rnap = d[d['shorthand']=='rnap']
    sigma = d[d['shorthand']=='all_sigma']
    ax.plot(rnap['growth_rate_hr'], rnap['n_complex'].values/sigma['n_complex'].values,
            'o', color=dataset_colors[g[0]], label=g[1], ms=4, markeredgewidth=0.5,
            markeredgecolor='k')

prot.viz.titlebox(ax, 'RNAP to σ-Factor Ratio', color='k', bgcolor=colors['pale_yellow'],
                  pad=0.05, boxsize=0.12, size=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/rnap_sigmafactor_ratio.pdf', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('total ribosome / total tRNA synthase', fontsize=6)

ax.set_ylim([0.5, 2.5])

for g, d in proteins.groupby(['dataset', 'dataset_name']):
    ribo = d[d['shorthand']=='ribosome']
    trna = d[d['shorthand']=='trna']
    ax.plot(trna['growth_rate_hr'], ribo['n_complex'].values/trna['n_complex'].values,
            'o', color=dataset_colors[g[0]], label=g[1], ms=4, markeredgewidth=0.5,
            markeredgecolor='k')

prot.viz.titlebox(ax, 'ribosome to tRNA synthase ratio', color='k', bgcolor=colors['pale_yellow'],
                  pad=0.05, boxsize=0.12, size=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/ribosome_trna_ratio.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('total ATP synthase /\ntotal ETC complexes', fontsize=6)

ax.set_ylim([0, 8])

for g, d in atp.groupby(['dataset', 'dataset_name']):
    _atp = d[d['shorthand']=='atp_synthase']
    proton = d[d['shorthand']=='proton_gradient']
    ax.plot(_atp['growth_rate_hr'], _atp['n_complex'].values/proton['n_complex'].values,
            'o', color=dataset_colors[g[0]], label=g[1], ms=4, markeredgewidth=0.5,
            markeredgecolor='k')

prot.viz.titlebox(ax, 'ATP synthase to proton pump ratio', color='k', bgcolor=colors['pale_yellow'],
                  pad=0.05, boxsize=0.12, size=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/atp_synthase_ETC_ratio.pdf', bbox_inches='tight')



# %%
