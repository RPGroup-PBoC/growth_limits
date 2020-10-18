#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import tqdm
# import upsetplot
from venn import venn
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}


# %%
# Load the compiled data sets
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# %%
#  Plot 1 ; total number of proteins quantified in each data set
fig, ax = plt.subplots(1,1, figsize = (4,3))

# for g, d in data.groupby('dataset'):

order = {'Schmidt et al. 2016' : 0,
         'Peebo et al. 2015' : 1,
         'Li et al. 2014' : 2,
         'Valgepea et al. 2013' : 3}

# creating the bar plot
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax.bar(order[g[1]], len(d['gene_name'].unique()),  color =dataset_colors[g[0]],
        width = 0.8)

plt.xticks([0,1,2,3])#, rotation='vertical')

ax.set_xticklabels(['Schmidt et al. 2016', 'Peebo et al. 2015', 'Li et al. 2014', 'Valgepea et al. 2013'], rotation=90)
# ax.set_xticks(rotation=75)
ax.set_ylabel('total number of\nproteins quantified')
plt.savefig('../../figures/intersections_venn_summed.pdf', bbox_inches='tight')

# %%
# Venn Diagram
fig = plt.figure()

datasets = {
    'Schmidt et al. 2016':  set(data[data['dataset']=='schmidt_2016']['gene_name'].unique()),
    'Peebo et al. 2015': set(data[data['dataset']=='peebo_2015']['gene_name'].unique()),
    'Li et al. 2014': set(data[data['dataset']=='li_2014']['gene_name'].unique()),
    'Valgepea et al. 2013' :  set(data[data['dataset']=='valgepea_2013']['gene_name'].unique())
}
venn(datasets, cmap = [colors['light_blue'], colors['green'], colors['purple'], colors['red']])

plt.savefig('../../figures/figA2_intersections_venn.pdf', bbox_inches='tight')


# %%
