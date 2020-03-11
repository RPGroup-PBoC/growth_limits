# %%
import pandas as pd
import tqdm
import glob

# Define the author names which did absolute measurments that are trustworth.
authors = ['schmidt2016', 'li2014', 'valgepea2013', 'peebo2015']

# Load inthe datasets and concatenate.
data = pd.concat([
    pd.read_csv(f'../../../data/{auth}_longform_annotated.csv') for auth in tqdm.tqdm(authors)],
    sort=False)
data.drop(columns=['reported_tot_per_cell', 'reported_fg_per_cell', 
                   'Unnamed: 0', 'reported_volume', 
                   'corrected_volume'], axis=1, inplace=True)

# Load the complexes to rename genes 
complexes = pd.read_csv('../../../data/ecocyc_raw_data/annotated_complexes.csv')

# Find synonyms in the data set and converge to give a single gene name 
no_synonyms = {}
for g, d in tqdm.tqdm(data.groupby('b_number')):
    gene = d['gene_name'].unique()[0]
    data.loc[data['b_number']==g, 'gene_name'] = gene 
    complexes.loc[complexes['b_number']==g, 'gene_name'] = gene 

# Rename the "not found" to "not assigned"
data.loc[data['cog_class']=='Not Found', 'cog_class'] = 'Not Assigned'
data.loc[data['cog_category']=='Not Found', 'cog_category'] = 'Not Assigned'

# Ignore the stationary state conditions. 
data = data[~data['condition'].str.contains('stationary_')]

# Save to disk
data.to_csv('../../../data/compiled_absolute_measurements.csv', index=False)

# Add counts to the molecular complexes
#%%
dfs = []
for g, d in tqdm.tqdm(complexes.groupby(['complex', 'b_number'])):
    # Find the number of subunits for each gene. 
    n_subunits = d['n_copies'].values[0]
    _df = data.loc[data['b_number']==g[1]].copy()
    _df['complex'] = g[0]
    _df['complex_annotation'] = d['annotation'].values[0]
    _df['n_subunits'] = n_subunits
    _df['n_units'] = _df['tot_per_cell'].values / n_subunits
    dfs.append(_df)

compiled_complexes = pd.concat(dfs, sort=False)
compiled_complexes.to_csv('../../../data/compiled_annotated_complexes.csv', index=False)
# %%
