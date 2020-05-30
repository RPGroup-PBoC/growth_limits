import numpy as np
import pandas as pd
from scipy import stats
import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Load the dataset with aboslute measurements
data = pd.read_csv('../../../data/compiled_datasets.csv')

# I want to plot the reported copy numbers from each dataset against
# Schmidt M9 glucose copy numbers.
# to keep things concise, only consider one growth condition from each dataset,
# selecting the condition that has a growth rate similar to Schmidt M9 glucose.
lambda_ref = data[data.dataset == 'schmidt_2016'][data.condition == 'glucose'].growth_rate_hr.unique()[0]

data['lambda_diff'] = abs(data['growth_rate_hr'] - lambda_ref)

data_trim = pd.DataFrame()
for d, df in data.groupby('dataset'):
    if d == 'schmidt_2016':
        df = df[df.condition == 'glucose']
    elif d == 'peebo_2015':
        df = df[df.condition == 'glucose_minimal']
        df = df[df.lambda_diff == df.lambda_diff.min()]
    else:
        df = df[df.lambda_diff == df.lambda_diff.min()]

    if d == 'li_2014':
        df['reported_tot_per_cell'] = df['tot_per_cell']

    data_trim = data_trim.append(df[['gene_name','reported_tot_per_cell', 'growth_rate_hr','dataset','dataset_name']])

data_trim = data_trim.reset_index()

d_schmidt = data_trim[data_trim.dataset == 'schmidt_2016'].copy()
d_li = data_trim[data_trim.dataset == 'li_2014'].copy()
d_ = [d_schmidt, d_li]

data_trim = data_trim[data_trim.dataset != 'schmidt_2016']
data_trim = data_trim[data_trim.dataset != 'li_2014']
N = len(data_trim.dataset.unique())

d_order = ['Valgepea et al. 2013', 'Peebo et al. 2015', 'Soufi et al. 2015',
        'Taniguichi et al. 2010']
d_dict = dict(zip(d_order, np.arange(N)))


## plotting!
fig, ax = plt.subplots(2, N, figsize=(2.5*N,6))

for d in d_order:
    data = data_trim[data_trim.dataset_name == d]

    ax[0, d_dict[d]].set_title(d,
        bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2})
    if d_dict[d] == 0:

        ax[0, d_dict[d]].set_ylabel('Schmidt et al. 2015\ncopies per cell',
            bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2})
        ax[1, d_dict[d]].set_ylabel('Li et al. 2014\ncopies per cell',
            bbox={'facecolor': '#EFCE9A', 'alpha': 0.5, 'pad': 2} )
    else:
        ax[0, d_dict[d]].get_yaxis().set_visible(False)
        ax[1, d_dict[d]].get_yaxis().set_visible(False)

    for i in [0,1]:
        ax[i, d_dict[d]].plot(np.logspace(0,7, 10),
                            np.logspace(0,7, 10),
                            ls = '-.',
                            alpha = 0.5,
                            c = 'k',
                            zorder = 0)

        _data = data.merge(d_[i], on='gene_name', suffixes=('', '_ref'))
        _data = _data.dropna()
        # we're going to scale log , so get rid of any item == zero
        _data = _data[_data.reported_tot_per_cell >0]
        _data = _data[_data.reported_tot_per_cell_ref >0]

        slope, intercept, r_value, _, _ = stats.linregress(np.log(_data.reported_tot_per_cell),
                                                        np.log(_data.reported_tot_per_cell_ref))

        ax[i, d_dict[d]].plot(np.exp(np.linspace(0.1,17, 100)),
                    np.exp(slope * np.linspace(0.1,17, 100) + intercept),
                    ls = '--',
                    alpha = 0.5,
                    zorder = 1,
                    label = '$r^2$: %.2f' % r_value**2)#+ np.round(str(r_value**2),decimals=2)) "%.2f" % a


        ax[i, d_dict[d]].scatter(_data.reported_tot_per_cell,
                              _data.reported_tot_per_cell_ref,
                              alpha = 0.5,
                              zorder = 10)

        ax[i, d_dict[d]].set_xlabel('copies per cell' )
        ax[i, d_dict[d]].legend()

for _ax in ax.ravel():
    if not _ax.lines:
        _ax.axis('off')
    else:
        _ax.set_xlim(1,1E6)
        _ax.set_ylim(1,1E6)
        _ax.set_xscale('log')
        _ax.set_yscale('log')



plt.tight_layout()
fig.savefig('../../../figures/dataset_correlations.pdf', bbox_inches='tight')
