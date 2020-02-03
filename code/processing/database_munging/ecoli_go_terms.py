"""
This script generates a dataframe which links individual E. coli genes with
their go terms. This generates a list that is highly degenerate with gene names
"""
#%%
import numpy as np 
import pandas as pd
import tqdm as tqdm

# Load the colicog list and GO term lookup table. 
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')
go_terms = pd.read_csv('../../../data/go_terms.csv')
go_terms.drop(labels='obsolete', axis=1, inplace=True)
go_terms.head()
# %%
# Instantiate a new dataframe.
df = pd.DataFrame([])
for g, d in tqdm.tqdm(colicogs.groupby('gene_name')):
    # Get the list of go terms
    terms = d['go_ids'].values[0].split('; ')
    for t in terms:
        _d = go_terms.loc[go_terms['id']==t]
        _d.rename(columns={'name':'go_term', 
                           'namespace':'go_namespace', 
                           'id':'go_number'}, inplace=True)
        _d['gene_name'] = g 
        df = df.append(_d)

# %%
df.to_csv('../../../data/escherichia_coli_go_terms.csv', index=False)

# %%
