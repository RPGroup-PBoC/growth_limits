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
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Load the original dataset with aboslute measurements
data_orig = pd.read_csv('../../../data/compiled_datasets.csv')

# Load the final dataset with aboslute measurements
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')

d_names = ['Schmidt et al. 2016', 'Li et al. 2014', 'Valgepea et al. 2013', 'Peebo et al. 2015']
d_dict = dict(zip(d_names, ['#7AA974', '#738FC1', '#D56C55', '#EAC264']))

## plotting!
fig, ax = plt.subplots(2, 2, figsize=(8,6))

for d, df in data_orig.groupby(['dataset_name']):
    if d not in d_names:
        continue
    for g, _df in df.groupby(['condition', 'growth_rate_hr']):
        if d == 'Li et al. 2014':
            mass = _df.fg_per_cell.sum()
        else:
            mass = _df.reported_fg_per_cell.sum()
        ax[0,0].scatter(_df.growth_rate_hr.unique()[0], mass,
                        color = d_dict[d], label = d)
        ax[0,0].set_xlabel('growth rate')
        ax[0,0].set_ylabel('total protein mass, fg')
        # volume prediction Si, F. et al. (2017), Current Biology, http://doi.org/10.1016/j.cub.2017.03.022
        # vol = 0.28 * np.exp(1.33  * _df.growth_rate_hr.unique()[0])
        vol = 0.27 * 2**(1.1  * _df.growth_rate_hr.unique()[0]/ np.log(2))
        conc = mass/ vol

        ax[0,1].scatter(_df.growth_rate_hr.unique()[0], conc,
                        color = d_dict[d], label = d)
        ax[0,1].set_xlabel('growth rate')
        ax[0,1].set_ylabel('total protein\nconcentration, fg/fL')


for d, df in data.groupby(['dataset_name']):
    if d not in d_names:
        continue
    for g, _df in df.groupby(['condition', 'growth_rate_hr']):

        mass = _df.fg_per_cell.sum()
        ax[1,0].scatter(_df.growth_rate_hr.unique()[0], mass,
                        color = d_dict[d], label = d)
        ax[1,0].set_xlabel('growth rate')
        ax[1,0].set_ylabel('total protein mass, fg')
        # volume prediction Si, F. et al. (2017), Current Biology, http://doi.org/10.1016/j.cub.2017.03.022
        # vol = 0.28 * np.exp(1.33  * _df.growth_rate_hr.unique()[0])
        vol = 0.27 * 2**(1.1  * _df.growth_rate_hr.unique()[0]/ np.log(2))
        conc = mass/ vol

        ax[1,1].scatter(_df.growth_rate_hr.unique()[0], conc,
                        color = d_dict[d], label = d)
        ax[1,1].set_xlabel('growth rate')
        ax[1,1].set_ylabel('total protein\nconcentration, fg/fL')

handles, labels = ax[1,0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax[1,0].legend(by_label.values(), by_label.keys(), loc = 'upper left', fontsize = 8)

# ax[1,0].set_ylim(0,1500)
ax[0,1].set_ylim(0,1100)
ax[1,1].set_ylim(0,1100)

ax[0,0].set_title('reported protein mass per cell',
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2})
ax[1,0].set_title('corrected protein mass per cell',
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2})
ax[0,1].set_title('reported protein concentration',
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2})
ax[1,1].set_title('corrected protein concentration',
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2})

plt.tight_layout()
fig.savefig('../../../figures/dataset_corrections.pdf', bbox_inches='tight')
