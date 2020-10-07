#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}

fig, ax = plt.subplots(1, 2, figsize=(6, 3))
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
# ax[1].set_yscale('log')
# ax[1].set_ylim([1E0, 1E3])
ax[0].set_ylabel('number of proteins per cell [million]', fontsize=6)
ax[1].set_ylabel('protein mass per cell [pg]', fontsize=6)
total_protein = data.groupby(['dataset', 'dataset_name', 
                                   'growth_rate_hr', 'condition'])['tot_per_cell'].sum().reset_index()
mass_protein = data.groupby(['dataset', 'dataset_name', 
                                   'growth_rate_hr', 'condition'])['fg_per_cell'].sum().reset_index()
                                   
for g, d in total_protein.groupby(['dataset', 'dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['tot_per_cell']/1E6, 'o', ms=4, color=dataset_colors[g[0]],
              label=g[1], markeredgewidth=0.5, markeredgecolor='k', alpha=0.75)
prot.viz.titlebox(ax[0], 'protein content', color='k', bgcolor=colors['pale_yellow'])
for g, d in mass_protein.groupby(['dataset', 'dataset_name']):
   ax[1].plot(d['growth_rate_hr'], d['fg_per_cell'] / 1E3, 'o', ms=4, color=dataset_colors[g[0]],
              label=g[1], markeredgewidth=0.5, markeredgecolor='k', alpha=0.75)
prot.viz.titlebox(ax[1], 'protein mass', color='k', bgcolor=colors['pale_yellow'])

for a in ax:
    a.legend(fontsize=6)
plt.savefig('../../figures/dataset_comparison.pdf', bbox_inches='tight')

# %%
sectors = data.groupby(['dataset', 'dataset_name', 'growth_rate_hr', 
                        'condition', 'cog_class'])['fg_per_cell'].sum().reset_index()

class_colors = {'information storage and processing': colors['blue'],
                'metabolism': colors['green'],
                'cellular processes and signaling': colors['red'],
                'Not Assigned': 'slategrey',
                'poorly characterized': 'slategrey'}

dataset_glyphs = {'li_2014': 'o', 'valgepea_2013': '^', 'peebo_2015':'v', 'schmidt_2016': 'd'}

fig, ax = plt.subplots(1, 2, figsize=(6, 3))
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('mass fraction per cell', fontsize=6)
ax[1].set_ylabel('metabolic mass fraction', fontsize=6)
ax[1].set_xlabel('information, storage, and processing\nmass fraction', fontsize=6)

for g, d in sectors.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
    tot_mass = d['fg_per_cell'].sum()
    for _g, _d in d.groupby(['cog_class']):
        ax[0].plot(_d['growth_rate_hr'], _d['fg_per_cell'].values/tot_mass, 
            linestyle='none', color=class_colors[_g], marker=dataset_glyphs[g[0]],
            label='__nolegend__', ms=4.5, markeredgewidth=0.5, markeredgecolor='k',
            alpha=0.75)
    _d = d[(d['cog_class']=='information storage and processing') |
           (d['cog_class']=='metabolism')]

    _d['fg_per_cell'] *= tot_mass**-1
    info = _d[_d['cog_class']=='information storage and processing']
    metab = _d[_d['cog_class']=='metabolism']
    ax[1].plot(info['fg_per_cell'], metab['fg_per_cell'], linestyle='none',
               color=dataset_colors[g[0]], marker='o', markeredgewidth=0.5,
               markeredgecolor='k', alpha=0.75, label='__nolegend__')

# Add legends. 
ax[0].plot([], [], 'ko', ms=4, label='Li et al. 2014')
ax[0].plot([], [], 'k^', ms=4, label='Valgepea et al. 2013')
ax[0].plot([], [], 'kv', ms=4, label='Peebo et al. 2015')
ax[0].plot([], [], 'kd', ms=4, label='Schmidt et al. 2016')
for s, v in class_colors.items():
    ax[0].plot([], [], '-', color=v, label=s.lower())
ax[1].plot([], [], 'o', markeredgewidth=0.5, markeredgecolor='k', 
            color=dataset_colors['li_2014'], ms=4, label='Li et al. 2014')
ax[1].plot([], [], 'o', markeredgewidth=0.5, markeredgecolor='k', 
            color=dataset_colors['valgepea_2013'], ms=4, label='Valgepea et al. 2013')
ax[1].plot([], [], 'o', markeredgewidth=0.5, markeredgecolor='k', 
            color=dataset_colors['peebo_2015'], ms=4, label='Peebo et al. 2015')
ax[1].plot([], [], 'o', markeredgewidth=0.5, markeredgecolor='k', 
            color=dataset_colors['schmidt_2016'], ms=4, label='Schmidt et al. 2016')

ax[0].legend(fontsize=6)
ax[1].legend(fontsize=6)
ax[0].set_ylim([0, 1.2])
plt.savefig('../../figures/sector_comparison.pdf', bbox_inches='tight')
# %%

# Look at individual genes. 
fig, ax = plt.subplots(3, 1, figsize=(5, 5), sharex=True)
prot.viz.titlebox(ax[0], 'elongation factor TufA', bgcolor=colors['pale_yellow'],
                size=6, color='k')
prot.viz.titlebox(ax[1], 'isocitrate lyase AceA', bgcolor=colors['pale_yellow'],
                size=6, color='k')
prot.viz.titlebox(ax[2], 'phosphofructokinase A PfkA', bgcolor=colors['pale_yellow'],
                size=6, color='k')



for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_ylabel('copy number per cell', fontsize=6)

ax[-1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

for g, d in data.groupby(['dataset', 'dataset_name']):
    tufA = d[d['gene_name']=='tufA']
    aceA = d[d['gene_name']=='aceA']
    pfkA = d[d['gene_name']=='pfkA']
    ax[0].plot(tufA['growth_rate_hr'], tufA['tot_per_cell'], 'o', ms=4.5, 
               color=dataset_colors[g[0]], markeredgewidth=0.5, 
               markeredgecolor='k', label=g[1])
    ax[1].plot(aceA['growth_rate_hr'], aceA['tot_per_cell'], 'o', ms=4.5, 
               color=dataset_colors[g[0]], markeredgewidth=0.5, 
               markeredgecolor='k', label=g[1])
    ax[2].plot(pfkA['growth_rate_hr'], pfkA['tot_per_cell'], 'o', ms=4.5, 
               color=dataset_colors[g[0]], markeredgewidth=0.5, 
               markeredgecolor='k', label=g[1])

for a in ax:
    a.set_yscale('log')
ax[0].set_ylim([1E4, 1E6])
ax[-1].set_ylim([1E2, 1E4])
ax[1].set_ylim([5E3, 1E6])
plt.savefig('../../figures/microscopic_agreement.pdf', bbox_inches='tight')
# %%
