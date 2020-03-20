#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import tqdm
import upsetplot
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# %%
# Load the compiled data sets
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')
set_df = pd.DataFrame([])
for g in tqdm.tqdm(data['gene_name'].unique()):
    d = data[data['gene_name']==g]['dataset'].unique()
    d = ', '.join(list(d))
    schmidt, li, peebo, valgepea = False, False, False, False
    if 'schmidt' in d:
        schmidt = True
    if 'li' in d:
        li = True 
    if 'valgepea' in d:
        valgepea = True 
    if 'peebo' in d:
        peebo = True
    set_df = set_df.append({'Schmidt et al. 2016': schmidt,
                            'Peebo et al. 2015': peebo,
                            'Li et al. 2014': li,
                            'Valgepea et al. 2013': valgepea,
                            'gene_name': 1},
                            ignore_index=True)

#%%
set_df.set_index([g for g in set_df.keys() if g != 'gene_name'], inplace=True)
#%% 
upsetplot.plot(set_df, sum_over='gene_name', show_counts=True, 
                sort_by='cardinality',
               facecolor=colors['blue'], element_size=45)
plt.savefig('../../figures/intersections.pdf', bbox_inches='tight')


# %%
