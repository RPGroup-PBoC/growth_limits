#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}

# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# calculate ribosomal fraction
L_R = 7459.0 # aa
ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']

ribo_frac_df = pd.DataFrame()
for c, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    mass_ribo = d[d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    frac_ribo = mass_ribo / d.fg_per_cell.sum()
    data_list = {'frac_ribo' : frac_ribo,
                'dataset' : c[0],
                'condition' : c[1],
                'growth_rate_hr' : c[2],
                'dataset_name' : d.dataset_name.values[0]}
    ribo_frac_df = ribo_frac_df.append(data_list,
                                        ignore_index = True)

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# Format and label the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlim([0, 1])
ax.set_ylim([0, (17.1 / L_R) * 3600 + 0.5])
ax.set_xlabel('ribosomal fraction ($\Phi_R$)', fontsize=6)
ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot the prediction.
frac = np.linspace(0,1.0,100)
gr = (17.1 * frac/ L_R) * 3600
ax.plot(frac, gr,  color='k', alpha=1.0, label='maximum growth rate',
                linestyle='-', lw=0.75)
ax.hlines((17.1 / L_R) * 3600, 0, np.max(gr), color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in ribo_frac_df.groupby(['dataset', 'dataset_name']):
    ax.plot(d['frac_ribo'], d['growth_rate_hr'], 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[1], ms=4)

ax.legend(ncol=1, fontsize=6)
plt.savefig('../../figures/ribosome_growth_limit.svg', bbox_inches='tight')
# %%
