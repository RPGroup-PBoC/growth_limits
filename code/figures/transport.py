#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')
data.head()

# %%
fig, ax = plt.subplots(2, 1, figsize=(4, 3))

_nitrogen = data[data['shorthand']=='nitrogen_tport']

ax[1].plot(_nitrogen['growth_rate_hr'], _nitrogen['n_complex'], 'o')
ax[1].set_ylim([0.1, 5000])
ax[1].set_yscale('log')

# %%
