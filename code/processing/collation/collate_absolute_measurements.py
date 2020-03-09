# %%
import pandas as pd
import glob

# Define the author names which did absolute measurments that are trustworth.
authors = ['schmidt2016', 'li2014', 'valgepea2013', 'peebo2015']

# Load inthe datasets and concatenate.
data = pd.concat([
    pd.read_csv(f'../../../data/{auth}_longform_annotated.csv') for auth in authors],
    sort=False)
data.drop(columns=['reported_tot_per_cell', 'reported_fg_per_cell', 
                   'annotation', 'Unnamed: 0', 'reported_volume', 
                   'corrected_volume'], axis=1, inplace=True)

# Find synonyms in the data set and converge to give a single gene name 
no_synonyms = {}
for g, d in data.groupby('b_number'):
    gene = d['gene_name'].unique()[0]
    data.loc[data['b_number']==g, 'gene_name'] = gene 

# Rename the "not found" to "not assigned"
data.loc[data['cog_class']=='Not Found', 'cog_class'] = 'Not Assigned'
data.loc[data['cog_category']=='Not Found', 'cog_category'] = 'Not Assigned'

# Ignore the stationary state conditions. 
data = data[~data['condition'].str.contains('stationary_')]

# Save to disk
data.to_csv('../../../data/compiled_absolute_measurements.csv', index=False)


# %%
