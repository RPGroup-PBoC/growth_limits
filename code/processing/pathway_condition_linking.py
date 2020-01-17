#%%
import numpy as np
import pandas as pd 

# Load the tidy schmidt data and uniprot pathway associations
condition_data = pd.read_csv('../../data/schmidt2016_longform.csv')
uniprot_data = pd.read_csv('../../data/uniprot_pathway_associations.csv')


# Set up the storage list for dataframes. 
dfs = []

# Iterate through each growth condition.
for g, d in condition_data.groupby('condition'):
    # Iterate through each pathway.
    for _g, _d in uniprot_data.groupby('pathway'):
        # Get only the genes associated with pathway in condition.
        pathway_elements = d[d['uniprot'].str.contains('|'.join(_d['uniprot'].values))]
       
# %%
