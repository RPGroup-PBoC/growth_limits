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
ppGpp = data[data['shorthand']=='ppGpp']

fig, ax = plt.subplots(1, 1, figsize=(3.5, 2))
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('RelA + SpoT', fontsize=6)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_ylim([0, 500])

for g, d in ppGpp.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', color=dataset_colors[g[0]],
            label=g[1], markeredgecolor='k', markeredgewidth=0.5, ms=4)

ax.legend(fontsize=6)
plt.savefig('../../figures/ppGpp_synthases.pdf')
# %%


# %%
