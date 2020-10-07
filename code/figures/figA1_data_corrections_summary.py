## plotting the total protein per cell and protein concentrations
# using the original data as reported, and our final compilated dataset.

## Here we are plotting the correlations across individual proteomic datasets

import numpy as np
import pandas as pd
from scipy import stats
import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import prot.viz
import prot.size
colors, palette = prot.viz.bokeh_theme()
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
         'peebo_2015':colors['green'], 'valgepea_2013':colors['red'],
         'soufi_2015':colors['yellow'], 'taniguichi_2010':colors['light_green']}
         # prot.viz.dataset_colors()
prot.viz.plotting_style()

# Load the original dataset with aboslute measurements
data_orig = pd.read_csv('../../data/compiled_datasets.csv')


# Load the final dataset with aboslute measurements
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

d_names = ['Schmidt et al. 2016', 'Li et al. 2014', 'Valgepea et al. 2013', 'Peebo et al. 2015',
        'Soufi et al. 2015', 'Taniguichi et al. 2010']
d_names2 = ['Schmidt et al. 2016', 'Li et al. 2014', 'Valgepea et al. 2013', 'Peebo et al. 2015']


d_dict = dict(zip(d_names, ['#7AA974', '#738FC1', '#D56C55', '#EAC264', '#905426', '#7AA974']))

## plotting!
fig, ax = plt.subplots(1, 2, figsize=(5,2.5))

for d, df in data_orig.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
    if d[1] not in d_names:
        continue
    if d[1] == 'Li et al. 2014':
        mass = df.fg_per_cell.sum()
    else:
        mass = df.reported_fg_per_cell.sum()
    ax[0].plot(df.growth_rate_hr.unique(), mass, 'o', ms=4, color=dataset_colors[d[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=d[1])


for d, df in data.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):#.groupby(['dataset_name']):

    if d[1] not in d_names2:
        continue

    mass = df.fg_per_cell.sum()

    ax[1].plot(df.growth_rate_hr.unique(), mass, 'o', ms=4, color=dataset_colors[d[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=d[1])



handles, labels = ax[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax[1].legend(by_label.values(), by_label.keys(), loc = 'upper left', fontsize = 6)

# ax[1,0].set_ylim(0,1500)
ax[0].set_ylim(0,850)
ax[1].set_ylim(0,850)

ax[0].set_title('reported protein mass per cell',
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2}, fontsize=8)
ax[1].set_title('final protein mass per cell',
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2}, fontsize=8)

for ax_ in ax:
    ax_.xaxis.set_tick_params(labelsize=5)
    ax_.yaxis.set_tick_params(labelsize=5)
    ax_.set_ylabel('total protein mass [fg]', fontsize=6)
    ax_.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

plt.tight_layout()
fig.savefig('../../figures/figA1_dataset_corrections.pdf', bbox_inches='tight')
