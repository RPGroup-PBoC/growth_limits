#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import prot.viz
import prot.size as size
colors, palette = prot.viz.bokeh_theme()

dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
prot.viz.plotting_style()


# %%
######################
# plot configuration #
######################
fig1, ax = plt.subplots(2, 1, figsize = (5,3))

# %%
######################
# plot 1  #
######################
# plot of avg copy num vs loc
# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# load in position information
regulonDB = pd.read_csv('../../data/regulonDB_raw/GeneProductSet.txt', delimiter = '	')

tss_map = dict(zip(regulonDB.b_number.values, regulonDB['Gene left end position in the genome'].values))
data['pos'] =data['b_number'].map(tss_map)

data_schmidt = data[data.dataset == 'schmidt_2016']
data_schmidt = data_schmidt.sort_values(by='growth_rate_hr', ascending = True)

###################
# Deal with double copies of EF-Tu
# Assume that they are equally distributed between
# gene copies. (right now, only tufA)
data_tufA = data_schmidt[data_schmidt.gene_name == 'tufA']
data_tufA['tot_per_cell'] = data_tufA['tot_per_cell']/2

data_tufB = data_schmidt[data_schmidt.gene_name == 'tufA']
data_tufB['gene_name'] = data_tufB['gene_name'].replace('tufA', 'tufB')
data_tufB['tot_per_cell'] = data_tufB['tot_per_cell']/2
data_tufB['pos'] = 4175944


data_schmidt = data_schmidt[data_schmidt.gene_name != 'tufA']
data_schmidt = data_schmidt.append(data_tufA)
data_schmidt = data_schmidt.append(data_tufB)

###################
###################

# data_schmidt_ = data_schmidt[data_schmidt.pos >= (4.6E6-0.25E6)]
data_schmidt_ = data_schmidt[data_schmidt.pos >= (4.6E6-1E6)]
data_schmidt_['pos'] = data_schmidt_['pos'] - 4.6E6
data_schmidt = data_schmidt.append(data_schmidt_)

# data_schmidt_ = data_schmidt[data_schmidt.pos <= (0.25E6)]
data_schmidt_ = data_schmidt[data_schmidt.pos <= (1E6)]
data_schmidt_['pos'] = data_schmidt_['pos'] + 4.6E6
data_schmidt = data_schmidt.append(data_schmidt_)

# color map stuff
# an array of parameters, each of our curves depend on a specific
# value of parameters
parameters = data_schmidt.growth_rate_hr.unique()

# norm is a class which, when called, can normalize data into the
norm = matplotlib.colors.Normalize(
    vmin=np.min(parameters),
    vmax=np.max(parameters))

c_m = matplotlib.cm.cividis_r
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)


# calculate avg
for c, d in data_schmidt.groupby('condition', sort= False):
    pos = np.linspace(0,4.6E6, 1000)
    avg_num = []
    for p in pos:
        d_ = d[d.pos <= p + 0.25E6]
        d_ = d_[d_.pos >= p - 0.25E6]
        avg_num = np.append(avg_num,d_.tot_per_cell.mean())

    # avg_num = avg_num - np.mean(avg_num)

    ax[0].plot((pos)/1E6, avg_num,
            color = s_m.to_rgba(d.growth_rate_hr.unique())[0],
            lw = 0.5)


# divider = make_axes_locatable(ax1)
# cax = divider.append_axes('right', size='5%', pad=0.05)
# cb = fig.colorbar(s_m, cax=cax, orientation='vertical');
#
# cb.ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
# cb.ax.tick_params(labelsize=5)
ax[0].set_xlabel('genomic position (Mb)', fontsize=6)
ax[0].set_ylabel('average protein copy\nnumber', fontsize=6)
ax[0].xaxis.set_tick_params(labelsize=5)
ax[0].yaxis.set_tick_params(labelsize=5)
ax[0].set_xlim(0,4.6)

# add in position info for oriC, ter
divider = make_axes_locatable(ax[0])
ax_genes = divider.append_axes('bottom', size='5%', pad=0.05)
ax_genes.set_xticklabels([])
ax_genes.set_yticklabels([])

ax_genes.set_xlim(0,4.6)
ax_genes.grid(False)
ax_genes.spines['bottom'].set_visible(False)
ax_genes.spines['right'].set_visible(False)
ax_genes.spines['top'].set_visible(False)
ax_genes.spines['left'].set_visible(False)
ax_genes.set_facecolor('white')
ax_genes.patch.set_alpha(0.0)

ax_genes.plot(np.array([1339769,1339769])/1E6, [-2,2], lw=2, color="#C24E9D", zorder = 10)
ax_genes.plot(np.array([1607181,1607181])/1E6, [-2,2], lw=2, color="#C24E9D", zorder = 10)
ax_genes.plot(np.array([3926090, 3926090])/1E6, [-2,2], lw=2, color="#19733A", zorder = 10)


####################
####################
# Plot 2 - remove ribosomal proteins and EF-Tu
ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO', 'tufA', 'tufB']

data_schmidt_null = data_schmidt[~data_schmidt.gene_name.isin(ribosome_genes)]
data_schmidt_null = data_schmidt_null.sort_values(by='growth_rate_hr', ascending = True)

# for c, d in data_schmidt_null.groupby('condition', sort= False):
#     pos = np.linspace(0,4.6E6, 1000)
#     avg_num = []
#     for p in pos:
#         d_ = d[d.pos <= p + 0.25E6]
#         d_ = d_[d_.pos >= p - 0.25E6]
#         avg_num = np.append(avg_num,d_.tot_per_cell.mean())
#
#     avg_num = avg_num - np.mean(avg_num)
#
#     ax[1].plot((pos)/1E6, avg_num,
#             color = s_m.to_rgba(d.growth_rate_hr.unique())[0],
#             lw = 0.5)

# divider = make_axes_locatable(ax2)
# cax = divider.append_axes('right', size='5%', pad=0.05)
# cb = fig.colorbar(s_m, cax=cax, orientation='vertical');
#
# cb.ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
# cb.ax.tick_params(labelsize=5)
ax[1].set_xlabel('genomic position (Mb)', fontsize=6)
ax[1].set_ylabel('average copy number\nrelative to mean', fontsize=6)
ax[1].xaxis.set_tick_params(labelsize=5)
ax[1].yaxis.set_tick_params(labelsize=5)
ax[1].set_xlim(0,4.6)

# add in position info for oriC, ter
divider = make_axes_locatable(ax[1])
ax_genes = divider.append_axes('bottom', size='5%', pad=0.05)
ax_genes.set_xticklabels([])
ax_genes.set_yticklabels([])

ax_genes.set_xlim(0,4.6)
ax_genes.grid(False)
ax_genes.spines['bottom'].set_visible(False)
ax_genes.spines['right'].set_visible(False)
ax_genes.spines['top'].set_visible(False)
ax_genes.spines['left'].set_visible(False)
ax_genes.set_facecolor('white')
ax_genes.patch.set_alpha(0.0)

ax_genes.plot(np.array([1339769,1339769])/1E6, [-2,2], lw=2, color="#C24E9D", zorder = 10)
ax_genes.plot(np.array([1607181,1607181])/1E6, [-2,2], lw=2, color="#C24E9D", zorder = 10)
ax_genes.plot(np.array([3926090, 3926090])/1E6, [-2,2], lw=2, color="#19733A", zorder = 10)

# fig1.savefig('../../figures/supplemental_boxcar1.pdf', bbox_inches = 'tight')

####################
####################
# Part 2 - effect of averaging window size - miller, 0 Mp

# %%
######################
# plot configuration #
######################
fig2, ax2 = plt.subplots(5, 1, figsize = (5,4))

# ax2[0].scatter((data_schmidt.pos.values)/1E6, (data_schmidt.tot.values),
#         color = s_m.to_rgba(d.growth_rate_hr.unique())[0],
#         ms = 1)

window_size = [0.05E6, 0.25E6, 0.5E6, 1.0E6, 2E6]
for i, ax in enumerate(ax2):
    # if i > 0:
    #     continue
    for c, d in data_schmidt_null.groupby('condition', sort= False):
        if i == 0:
            pos = np.linspace(0,4.6E6, 10000)
        else:
            pos = np.linspace(0,4.6E6, 1000)
        avg_num = []
        pos_ = []
        for p in pos:
            d_ = d[d.pos <= p + window_size[i]/2]
            d_ = d_[d_.pos >= p - window_size[i]/2]
            if d_.empty:
                continue
            pos_ = np.append(pos_,p)
            avg_num = np.append(avg_num,d_.tot_per_cell.mean())
        
        avg_num = avg_num - np.mean(avg_num)

        ax.plot((pos_)/1E6, avg_num,
                color = s_m.to_rgba(d.growth_rate_hr.unique())[0],
                lw = 0.5)

        ax.set_xlabel('genomic position (Mb)', fontsize=6)
        # ax.set_ylabel('average copy number\nrelative to mean', fontsize=6)
        ax.xaxis.set_tick_params(labelsize=5)
        ax.yaxis.set_tick_params(labelsize=5)
        ax.set_xlim(0,4.6)

fig2.savefig('../../figures/supplemental_boxcar4.pdf', bbox_inches = 'tight')
