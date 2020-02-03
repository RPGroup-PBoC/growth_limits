"""
This script takes a dataframe of absolute protein measurements and breaks
each gene down by its annotated go term. This results in a large dataframe that 
is highly redundant, but will be useful for generating plots of molecular 
complex abundance as a function of growth rate. 

TODO: Figure out how to properly deal with quaternary structure and subunit
abundance
"""
#%%
import numpy as np
import pandas as pd 
import tqdm as tqdm

# Load the compiled dataset of absolute measurements
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

# Load the go-ter look-up-table
go_term = pd.read_csv('../../data/escherichia_coli_go_terms.csv')

# Merge the two on the gene namees. 
complete = data.merge(go_term, on='gene_name')

# Save to disk
if 'Unnamed: 0' in complete.keys():
    complete.drop(labels='Unnamed: 0', axis=1, inplace=True)
if 'growth_rate' in complete.keys():
    complete.drop(labels='growth_rate', axis=1, inplace=True)
complete.to_csv('../../data/compiled_absolute_go_terms.csv', index=False)
