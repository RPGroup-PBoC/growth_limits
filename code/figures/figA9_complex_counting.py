#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors = prot.viz.plotting_style()

# Load the complex data and find the largest 6
data = pd.read_csv('../../data/compiled_annotated_complexes.csv')
data = data[(data['dataset']=='schmidt_2016') & (data['condition']=='glucose')]
nums = data.groupby(['complex', 'complex_annotation'])['gene_name'].count().reset_index()
nums.sort_values(by='gene_name', inplace=True)

# Define the complexes to show
cplxs = ['CPLX0-3964', 'NADH-DHI-CPLX', 'ATPSYN-CPLX',
         'CPLX0-7992', 'SEC-SECRETION-CPLX', 'CPLX0-8301']
# %%
fig, ax = plt.subplots(3, 2, figsize=(7.5, 6))
ax = ax.ravel()
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_ylabel(r'$\frac{\mathrm{observed\,copy\,number}}{\mathrm{annotated\,number\,in\,complex}}$', fontsize=8)

for i, cplx in enumerate(cplxs):
    prot_cplx = data[data['complex']==cplx]
    prot.viz.titlebox(ax[i], prot_cplx['complex_annotation'].values[0],
                     bgcolor=colors['pale_yellow'], color='k')

    # Compute the mean complex number.
    mean_num = prot_cplx['n_units'].mean()
    ax[i].hlines(mean_num, -1, len(prot_cplx), lw=3, color='grey', label='arithmetic mean',
                alpha=0.3)
    ax[i].plot(np.arange(0, len(prot_cplx), 1), prot_cplx['n_units'], 'o',
                    markeredgewidth=0.5, markeredgecolor='k',
                    color=colors['dark_green'], label='__nolegend__', ms=4)

    ax[i].set_xlim([-.5, len(prot_cplx) - .5])
    ax[i].set_xticks(np.arange(0, len(prot_cplx), 1))
    labels = [g[0].upper() + g[1:] for g in prot_cplx['gene_name'].values]
    ax[i].set_xticklabels(labels)
    ax[i].set_yscale('log')
    ax[i].set_ylim([10**(np.log10(prot_cplx['n_units'].min()) - 1),
                   10**(np.log10(prot_cplx['n_units'].max()) + 1)])
    if i == 0:
        ax[i].xaxis.set_tick_params(labelsize=3)
        ax[i].set_xticklabels(labels, rotation=90)
ax[0].legend(loc='lower right', fontsize=6)
plt.savefig('../../figures/figA9_subunit_counting.pdf', bbox_inches='tight')
