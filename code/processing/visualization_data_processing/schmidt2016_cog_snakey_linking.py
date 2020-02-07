"""
This script generates a csv amenable to visualization using a snakey diagram.
The nodes are definde by the cog category, cog class, and gene name with the
width of the edges related to the mass fraction from the parent node to child. 
"""
#%%
import numpy as np
import pandas as pd

# Load the schimdt long-form dataset. 
data = pd.read_csv('../../../data/schmidt2016_longform_annotated.csv')

# Make a unique snakey graph for each condition
for g, d in data.groupby(['condition']):

    # Compute the total proteome mass *per condition* avoiding duplicate entries
    total_mass = d[~d.duplicated(subset='gene_name')]['fg_per_cell'].sum()

    # Compute the fractional mass for each  

# %%


# %%
