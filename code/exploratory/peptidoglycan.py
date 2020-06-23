#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

# Load the data and restrict. 
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
tpds = data[data['shorthand']=='transpeptidases']

fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.set_yscale('log')
ax.set_ylim([1E0, 1E4])


theta_pg = 0.025 * 0.3
rho = constants['density']['value']
cell_mass = constants['cell_mass']['value']
m_pg = 1E3 * 1E12 / 6E23 
kcat = 2 
n_tpds = (theta_pg * cell_mass /2 ) / (m_pg * kcat * constants['t_double']['value'])
ax.plot(constants['growth_rate']['value'], n_tpds, '-', color='grey', alpha=0.5,
        label='cell size dependence')
for g, d in tpds.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, label=g[1], color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k')

ax.legend(fontsize=6)
# %%
