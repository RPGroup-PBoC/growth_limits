#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import scipy
import seaborn as sns
import tqdm
import prot.viz
import prot.size as size
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()

# Load the datasets.
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')
regulonDB = pd.read_csv('../../data/regulonDB_raw/GeneProductSet.txt', delimiter = '	')

# Map the transcription start stie to position from regulon DB.
tss_map = dict(zip(regulonDB.b_number.values, regulonDB['Gene left end position in the genome'].values))
data['pos'] =data['b_number'].map(tss_map)

data_schmidt = data[data.dataset == 'schmidt_2016']
data_schmidt = data_schmidt.sort_values(by='growth_rate_hr', ascending = True)

# Deal with double copies of EF-Tu
# Assume that they are equally distributed between gene copies. (right now, only tufA)
data_tufA = data_schmidt[data_schmidt.gene_name == 'tufA']
data_tufA['tot_per_cell'] = data_tufA['tot_per_cell']/2

data_tufB = data_schmidt[data_schmidt.gene_name == 'tufA']
data_tufB['gene_name'] = data_tufB['gene_name'].replace('tufA', 'tufB')
data_tufB['tot_per_cell'] = data_tufB['tot_per_cell']/2
data_tufB['pos'] = 4175944

data_schmidt = data_schmidt[data_schmidt.gene_name != 'tufA']
data_schmidt = data_schmidt.append(data_tufA)
data_schmidt = data_schmidt.append(data_tufB)

#  Shift the genomic positions to be centered at oriC
oriC_loc = 3925744
data_schmidt['shifted_pos'] = data_schmidt['pos'].values - oriC_loc
max_pos_noshift = data_schmidt['pos'].max() + 1
max_pos = data_schmidt['shifted_pos'].max()
data_schmidt.loc[data_schmidt['shifted_pos'] < 0, 'shifted_pos'] += max_pos_noshift
data_schmidt.dropna(inplace=True)

# Create a mapping between shifted position and radians
pos_to_rad = {p:v for p, v in zip(np.sort(data_schmidt['shifted_pos'].unique()), np.linspace(0, 2 * np.pi, len(data_schmidt['shifted_pos'].unique())))}

# Define the color palette to be used
palette = sns.color_palette('viridis', n_colors=len(data['growth_rate_hr'].unique()))
color_mapping = {k:v for k, v in zip(data_schmidt['growth_rate_hr'].unique(), palette)}




#%%

# %%
######################
# plot configuration #
######################
fig = plt.figure(constrained_layout=True)
# widths = [6, 2.5, 2.5, 5]
widths = [6, 5, 5]
heights = [1.5, 1, 0.75, 1.25, 0.75, 2]
spec = fig.add_gridspec(ncols=3, nrows=6, width_ratios=widths,
                          height_ratios=heights)

# # plot of t_cyc vs. tau
# ax1 = fig.add_subplot(spec[0, 1])
# # plot of t_cyc vs. tau
# ax2 = fig.add_subplot(spec[0, 2])
# # plot of RNA/protein vs. num ori
# ax3 = fig.add_subplot(spec[1:4, 1])
# # plot of # ribosomes vs. num ori
# ax4 = fig.add_subplot(spec[1:4, 2])
# # plot of avg copy num vs loc
ax5 = fig.add_subplot(spec[3, 0])


# plot of elongation rate vs. active fraction
ax6 = fig.add_subplot(spec[4:, 1])
ax6_ = fig.add_subplot(spec[4:, 0])
ax6_.axis('off')

# # plot of ribosomal fraction vs lambda
# ax7 = fig.add_subplot(spec[4:, 2])

######################
# plot of avg copy num vs loc
# # Load the compiled data
# Deal with edge positions
# re-copy values with positions near
# 0 and 4.6E6 so that we can treat as circular chromosome.
# 0 bp side
data_schmidt_ = data_schmidt[data_schmidt.pos >= (4.6E6-0.5E6)]
data_schmidt_['pos'] = data_schmidt_['pos'] - 4.6E6
data_schmidt = data_schmidt.append(data_schmidt_)
# 4.6E6 bp side
data_schmidt_ = data_schmidt[data_schmidt.pos <= (0.5E6)]
data_schmidt_['pos'] = data_schmidt_['pos'] + 4.6E6
data_schmidt = data_schmidt.append(data_schmidt_)

# colormap stuff
colors_viridis = plt.cm.cividis(np.linspace(0,1,len(data_schmidt.condition.unique())))
colordic = dict(zip(data_schmidt.condition.unique(), colors_viridis))
parameters = data_schmidt.growth_rate_hr.unique()

#  normalize data for colormap range from 0 to 2 hr-1
norm = matplotlib.colors.Normalize(
    vmin=np.min(parameters),
    vmax=np.max(parameters))

c_m = matplotlib.cm.cividis_r
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)

###################
# Now calculate average values for plotting!
for c, d in data_schmidt.groupby('condition', sort= False):
    pos = np.linspace(0,(4.6E6), 500)
    avg_num = []
    for p in pos:
        # d_ = d[d.pos <= p + 50000]
        # d_ = d_[d_.pos >= p - 50000]

        # create weighting
        positions = d.pos.values
        weights = scipy.stats.norm(p, 0.1E6).pdf(positions)

        averages = np.sum(d.tot_per_cell.values * weights)

        avg_num = np.append(avg_num, averages)

    avg_num = avg_num - np.mean(avg_num)

    ax5.plot((pos)/1E6, avg_num,
            color = s_m.to_rgba(d.growth_rate_hr.unique())[0],
            lw = 0.5)


divider = make_axes_locatable(ax5)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(s_m, cax=cax, orientation='vertical');

cb.ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
cb.ax.tick_params(labelsize=5)
ax5.set_xlabel('genomic position (Mb)', fontsize=6)
ax5.set_ylabel('average protein copy\nnumber relative to mean', fontsize=6)
ax5.xaxis.set_tick_params(labelsize=5)
ax5.yaxis.set_tick_params(labelsize=5)
# ax5.set_xlim(0,4.6E6)
ax5.set_xlim(0,4.6)

# add in position info for rRNA, r-protein, oriC, ter
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

rrna_operons = np.array([4035239, 4166367, 3941516, 3429047, 4207863, 2731448, 223593])
for p in rrna_operons:
    ax_genes.plot(np.array([p,p])/1E6, [-2,2], lw=1, color="#E55E68")

pos_schmidt = data[data.dataset == 'schmidt_2016']
data['pos'] =data['b_number'].map(tss_map)
ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']
# ribosome_genes = ['tufA']

pos_schmidt = pos_schmidt[pos_schmidt.gene_name.isin(ribosome_genes)].pos.unique()
for p in pos_schmidt:
    ax_genes.plot(np.array([p,p])/1E6, [-2,2], lw=0.5, color="#F47530")

fig.savefig('../../figures/gaussian_avg_copynumber.pdf', bbox_inches = 'tight')
