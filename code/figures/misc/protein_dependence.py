#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')
trna = data[data['shorthand']=='trna']
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}


# %%
rate_range = np.linspace(0.01, 2, 200)
t_double = 3600 * np.log(2) / rate_range
rho = 1.1 # in pg/um^3
vol = prot.size.lambda2size(rate_range)
phi_p = 0.15
m_aa = 110 / (6.022E23 * 1E-12) # in pg
k_cat = 20
N_synthases = (rho * vol * phi_p) / (m_aa * t_double * k_cat)

fig, ax = plt.subplots(1,1)
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.set_ylim([1E3, 1E6])
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('number of tRNA sythetases')

ax.plot(rate_range, N_synthases, 'k--', label='prediction')
for g, d in trna.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', color=dataset_colors[g[0]],
            ms=4.5, markeredgewidth=0.75, markeredgecolor='k', label=g[1])
ax.legend()
plt.savefig('../../figures/tRNA_synthase_dependence.pdf')
# %%
