#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()


data = pd.read_csv('../../data/compiled_estimate_categories.csv')
dntp = data[data['shorthand']=='dntp']
dnap = data[data['shorthand']=='dnap']
sig = data[data['shorthand']=='sigma70']
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}



# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])

# Define the estimate
L_genes = 4500
L_spacing = 80
t_div = (2.5E6 / 600) / (60 * 60)

rate_range = np.linspace(0, 2, 200)
prediction = (L_genes / L_spacing) * 7 * 2**(t_div * rate_range)

ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in sig.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('RpoD ($\sigma_{70}$) per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/sigmafactor_trend.pdf', bbox_inches='tight')

# %%
# predicted_dnap
prediction  = (1E7 * 2**(t_div * rate_range)) /(10 * 60 * 60 * rate_range**-1)
prediction[prediction < 1] = 1
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
# ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])


ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in dntp.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribonucleoreductases per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/rnr_trend.pdf', bbox_inches='tight')



# %%
prediction  = (5E6 * 2**(t_div * rate_range)) / (60 * 60 * 600 * rate_range**-1)
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
# ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])


ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in dnap.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribonucleoreductases per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/rnr_trend.pdf', bbox_inches='tight')



# %%
