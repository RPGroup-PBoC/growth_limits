#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/compiled_annotated_complexes.csv', comment='#')
data = data[data['growth_rate_hr']==0.5]
ribo = data[data['complex']=='CPLX0-3964']
dnap = data[data['complex']=='CPLX0-3803']
atp = data[data['complex']=='ATPSYN-CPLX']

dnap_locs = {c:i for i, c in enumerate(dnap['gene_name'].unique())}
ribo_locs = {c:i for i, c in enumerate(ribo['gene_name'].unique())}
atp_locs = {c:i for i, c in enumerate(atp['gene_name'].unique())}

# Compute the mean  and difference from mean
for df, locs in zip([dnap, ribo, atp], [dnap_locs,  ribo_locs, atp_locs]):
    mean_no = df['n_units'].mean()
    df['diff'] = df['n_units'] - mean_no
    for sub, pos in locs.items():
        df.loc[df['gene_name']==sub, 'pos'] = pos

#%%
fig, ax = plt.subplots(3, 1, figsize=(6, 6))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=6)
    a.hlines(0, -1, 55, 'k', linestyle='--')
    a.set_ylabel('functional subunit abundance\ndifference from mean', fontsize=6)

# Format the axes and add labels
ax[0].set_xlabel('ribosomal subunit', fontsize=6)
ax[1].set_xlabel('DNA polymerase III subunit',  fontsize=6)
ax[2].set_xlabel('F$_1$-F$_0$ ATP synthase subunit',  fontsize=6)
ax[0].set_xlim([-1, 55])
ax[1].set_xlim([-1, 7])
ax[2].set_xlim([-1, 8])

# Add titles
prot.viz.titlebox(ax[0], 'mature ribosome complex', fontsize=6, color='k',
                  bgcolor=colors['pale_yellow'], boxsize=0.12)
prot.viz.titlebox(ax[1], 'DNA polymerase III core enzyme', fontsize=6, color='k',
                  bgcolor=colors['pale_yellow'], boxsize=0.12)
prot.viz.titlebox(ax[2], 'F$_1$-F$_0$ ATP Synthase', fontsize=6, color='k',
                  bgcolor=colors['pale_yellow'], boxsize=0.12)


for i, l in enumerate([ribo_locs, dnap_locs, atp_locs]):
    _ = ax[i].set_xticks(list(l.values()))
    if i == 0:
        rot = 90
    else:
        rot = 0
    _ = ax[i].set_xticklabels(list(l.keys()), rotation=rot)


ax[0].plot(ribo['pos'], ribo['diff'], 'o', color=colors['dark_green'],
        alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', ms=4)
ax[1].plot(dnap['pos'], dnap['diff'], 'o', color=colors['dark_green'],
        alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', ms=4)
ax[2].plot(atp['pos'], atp['diff'], 'o', color=colors['dark_green'],
        alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', ms=4)
plt.tight_layout()
plt.savefig('../../figures/subunit_differences.pdf', bbox_inches='tight')
# %%
