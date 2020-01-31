#%%
import numpy as np
import pandas as pd
import tqdm

# Load the complete gene-linked dataset
data = pd.read_csv('../../data/go_biological_processes.csv')

# %%
# Make a dictionary of all possible combinations
processes = data['go_biological_process'].unique()

connection_dict = {f'{p1}:::{p2}':0 for p1 in processes for p2 in processes if p1 != p2} 


# %%
# Iterate through each gene and populate the connections
for g, d in tqdm.tqdm(data.groupby('gene')):
    # Find the processes. 
    member_of = d['go_biological_process'].values
    for p1 in member_of:
        for p2 in member_of:
            if p1 != p2:
                connection_dict[f'{p1}:::{p2}'] += 1
