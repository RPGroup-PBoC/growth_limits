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
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
rho = constants['density']['value']
cell_mass = constants['cell_mass']['value']

m_pg = 1E3 / 6E11 # in pg
kcat = 2
l = prot.size.lambda2length(growth_rate)
w = prot.size.lambda2width(growth_rate)
v = prot.size.lambda2size(growth_rate)
w_pg = 0.005 # in µm

# Determine the density. 
V_std = prot.size.rod_SA(2, 1, 1) * w_pg
mass_pg = 0.025 * 0.3 # in pg
rho_pg = mass_pg / V_std
V_pg = prot.size.rod_SA(l, w, v) * w_pg

n_tpds = (0.3 * V_pg * rho_pg) / (m_pg * kcat * t_double)
# n_tpds = (theta_pg * cell_mass) / (m_pg * kcat * constants['t_double']['value'])
ax.plot(constants['growth_rate']['value'], n_tpds, '-', color='grey', alpha=0.5,
        label='cell size dependence')
for g, d in tpds.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, label=g[1], color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k')

ax.legend(fontsize=6)
plt.savefig('./pg.pdf')
# %%
