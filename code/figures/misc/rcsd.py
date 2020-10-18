#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import prot.viz 
import prot.size
_ = prot.viz.plotting_style()


colors = prot.viz.dataset_colors()

data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
data.head()
# %%
rsd = data[data['gene_name']=='rsd']



fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.set_xlabel('growth rate [hr$^{1}$]')
ax.set_ylabel('Rsd copy number')
for g, d in rsd.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['tot_per_cell'], 'o', color=colors[g[0]],
            label=g[1], markeredgecolor='k', markeredgewidth=0.25,
            ms=8)

ax.legend(fontsize=8)
plt.savefig('./rsd_copy_number.pdf')



# %%  Compute the concentration

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.set_xlabel('growth rate [hr$^{1}$]')
ax.set_ylabel('Rsd concentration [µM]')
for g, d in rsd.groupby(['dataset', 'dataset_name']):
    vol = prot.size.lambda2size(d['growth_rate_hr']) * 1E-15
    conc = (d['tot_per_cell'].values / 6.022E23) / vol 
    ax.plot(d['growth_rate_hr'], conc * 1E6, 'o', color=colors[g[0]],
            label=g[1], markeredgecolor='k', markeredgewidth=0.25,
            ms=8)

ax.legend(fontsize=8)
plt.savefig('./rsd_concentration.pdf')


# %%
