#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import seaborn as sns
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

sns.set_palette('tab20b', n_colors=30)
# Load the data 
data = pd.read_csv('../../data/compiled_annotated_complexes.csv')

# %%
sugar_tporters = data[data['go_terms'].str.contains('GO:0008643')]

# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.set_yscale('log')
ax.set_ylabel('number of complexes')
ax.set_xlabel('growth rate [hr$^{-1}$]')
prot.viz.titlebox(ax, 'carbohydrate transporters', color='k', bgcolor=colors['pale_yellow'])

for g, d in sugar_tporters.groupby(['complex_annotation']):
    g = g.replace('<I>', '$')
    g = g.replace('<i>', '$')
    g = g.replace('</i>', '$')
    g = g.replace('</I>', '$')
    g = g.replace('&beta;', 'β')
    _d = d.groupby(['condition', 'growth_rate_hr'])['n_units'].mean().reset_index()
    ax.plot(_d['growth_rate_hr'], _d['n_units'], 'o', label=g, alpha=0.2, 
            markeredgewidth=0.25, markeredgecolor='k')
_sugar_tporters = sugar_tporters.groupby(['growth_rate_hr', 'dataset', 
                            'complex_annotation'])['n_units'].mean().reset_index()
_sugar_tporters = _sugar_tporters.groupby(['growth_rate_hr'])['n_units'].sum().reset_index()
ax.plot(_sugar_tporters['growth_rate_hr'], _sugar_tporters['n_units'], 'o', 
        color=colors['red'], markeredgecolor='k', markeredgewidth=0.5, 
        label='sum total')
ax.legend(bbox_to_anchor=(1.02, 1.02), ncol=1, fontsize=6)
plt.savefig('../../figures/all_sugar_transporters.pdf', bbox_inches='tight')
# %%


# %%
