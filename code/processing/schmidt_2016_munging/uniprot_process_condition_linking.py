#%%
import numpy as np
import pandas as pd 

# Load the tidy schmidt data and uniprot pathway associations
condition_data = pd.read_csv('../../data/schmidt2016_longform.csv')
uniprot_data = pd.read_csv('../../data/uniprot_biological_processes.csv')

# Set up the storage list for dataframes. 
dfs = []

# Iterate through each growth condition.
for g, d in condition_data.groupby('condition'):

    # Iterate through each pathway.
    for _g, _d in uniprot_data.groupby('process'):
        # Get only the genes associated with pathway in condition.
        process_members = d[d['uniprot'].str.contains('|'.join(_d['uniprot'].values))]

        # Add identifier of the biological process
        process_members['uniprot_bio_process'] = _g
        dfs.append(process_members)
linked_processes = pd.concat(dfs, sort=False)
linked_processes.to_csv('../../data/schmidt2016_genes_processes.csv', index=False)
# %%
# Compute the fraction information.
process_df = pd.DataFrame([])
for g, d in linked_processes.groupby(['condition', 'uniprot_bio_process']):
    # Get the total mass and size of the proteome in this condition. 
    proteome_mass = condition_data[
                    condition_data['condition']==g[0]]['fg_per_cell'].sum()
    proteome_size = condition_data[
                    condition_data['condition']==g[0]]['tot_per_cell'].sum()

    # Compute the fractioning
    mass_frac = d['fg_per_cell'].sum() / proteome_mass
    count_frac = d['tot_per_cell'].sum() / proteome_size

    # Append information to the dataframe. 
    process_dict = {'condition':g[0], 
                    'uniprot_bio_process':g[1],
                    'frac_mass': mass_frac,
                    'frac_count': count_frac}
    process_df = process_df.append(process_dict, ignore_index=True)

process_df.to_csv('../../data/proteome_process_sectoring.csv', index=False)

# %%
