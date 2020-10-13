#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
import prot.viz
import prot.size as size
import scipy.stats
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
# pos_to_rad = {p:v for p, v in
# zip(np.sort(data_schmidt['shifted_pos'].unique()), np.linspace(0, 2 * np.pi,
# len(data_schmidt['shifted_pos'].unique())))}
_positions = np.arange(0, 4.6E6, 1)
pos_to_rad = {p:v for p, v in zip(_positions, np.linspace(0, 2 * np.pi, len(_positions)))}


# Define the color palette to be used
palette = sns.color_palette('viridis', n_colors=len(data['growth_rate_hr'].unique()))
color_mapping = {k:v for k, v in zip(np.sort(data_schmidt['growth_rate_hr'].unique()), palette)}

#
ribo_prots = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']


# Set up the figure canvas,
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_zero_location('N')
ax.set_ylim([-60, 80])
# ax.set_yticks([])
ax.set_rlabel_position(180)
# ax.set_yscale('log')

# Play with one condition
cond = data_schmidt[data_schmidt['condition']=='42C']
cond.sort_values(by='shifted_pos', inplace=True)
window =  2E4
position = np.sort(data_schmidt['shifted_pos'].unique())
z = 100
iter =0
for g, d in tqdm.tqdm(data_schmidt.groupby(['growth_rate_hr']), desc='processing growth rates'):
    cond = d
    means = []
    pos = []
    for p in _positions[::500]:
        weights = scipy.stats.norm(p, window).pdf(d['shifted_pos'].values)
        means.append(np.sum(weights * d['tot_per_cell'].values))
        pos.append(pos_to_rad[p])

    shifted = list(np.array(means) - np.median(means))
    shifted.append(shifted[0])
    pos.append(pos[0])
    ax.plot(pos, shifted, '-', lw=1, alpha=0.8,  zorder=z, color=palette[iter])
    z -= 1
    iter += 1

ribosomes = data_schmidt[data_schmidt['gene_name'].isin(['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO'])]

for g, d in ribosomes.groupby(['pos', 'shifted_pos']):
    ax.plot([pos_to_rad[g[1]], pos_to_rad[g[1]]], [-15, -5], '-', lw=1, color=colors['red'])

plt.savefig('../../figures/fig11b_polar_chromosome.svg')
#%%
