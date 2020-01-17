#%%
import pandas as pd 
import numpy as np

#Load the raw schmidt data to generate a dictionary of uniprot accession to name
data = pd.read_csv('../../data/schmidt2016_longform.csv')

# Load the raw data 
uniprot_data = pd.read_csv('../../data/schmidt2016_uniprot.tab', delimiter='\t')

# Load the keyword id list. 
keywords = pd.read_csv('../../data/uniprot_keywords.csv')

# Restrict to the biological processes
processes = keywords[keywords['category']=='Biological process']
#%%
# Load the keyword ids.
gene_name_dict = {}
cog_class_dict = {}
cog_letter_dict = {}
cog_desc_dict = {}
gene_desc_dict = {}
for g, d in data.groupby(['uniprot']):
    gene_name_dict[g] = d['gene'].values[0]
    cog_desc_dict[g] = d['cog_desc'].values[0]
    cog_letter_dict[g] = d['cog_class_letter'].values[0]
    cog_class_dict[g] = d['cog_class'].values[0]
    gene_desc_dict[g] = d['desc'].values[0]

# %%
#%%
# Rename some columns to be more meaningful.
uniprot_data.rename(columns={'Interacts with': 'interactions', 
                     'Gene ontology (biological process)':'go_process',
                     'Gene ontology (molecular function)':'go_function',
                     'Cofactor':'cofactor'}, inplace=True)

# Generate an empty dataframes for biological process and cofactors
process_df = pd.DataFrame([])

# Group by protein and split each entry as appropriate. 
for g, d in uniprot_data.groupby(['search']):
    # Deal with a mistake in the uniprot accession info
    if g == 'P00452,P00452-2':         
        g = 'P00452'
    if g == 'P02919,P02919-2':
        g = 'P02919'

    # Get the list of keywords
    keys = d['Keyword ID'].values[0]
    if str(keys) != 'nan':
        keyIDs = d['Keyword ID'].values[0].split('; ')
    else:
        keyIDS = []
    # Iterate through each keyword ID, determine if it's a process, and store.  
    for key in keyIDs:
        val = keywords[keywords['accession_number']==key]
        if val['category'].values[0].lower() == 'biological process':
            # Get the process, cog names, and other inportant info.
            process_dict = {
                'process':val['id'].values[0],
                'uniprot':g,
                'gene': gene_name_dict[g],
                'desc': gene_desc_dict[g],
                'cog_desc':  cog_desc_dict[g],
                'cog_class': cog_class_dict[g],
                'cog_class_letter': cog_letter_dict[g]
                }
            process_df = process_df.append(process_dict, ignore_index=True)
process_df.to_csv('../../data/uniprot_biological_processes.csv', index=False)