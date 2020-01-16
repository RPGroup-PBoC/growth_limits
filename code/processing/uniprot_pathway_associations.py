#%%
import pandas as pd 
import numpy as np

#Load the raw schmidt data to generate a dictionary of uniprot accession to name
data = pd.read_csv('../../data/schmidt2016_dataset.csv')
gene_name_dict = {}
cog_class_dict = {}
cog_letter_dict = {}
cog_desc_dict = {}
gene_desc_dict = {}
for g, d in data.groupby(['Uniprot Accession']):
    gene_name_dict[g] = d['Gene'].values[0]
    cog_class_dict[g] = d['Annotated functional COG group (description)'].values[0]
    cog_letter_dict[g] = d['Annotated functional COG groups (letter)'].values[0]
    cog_desc_dict[g] = d['Annotated functional COG class'].values[0]
    gene_desc_dict[g] = d['Description'].values[0].split('OS=')[0][:-1]


# %%
# Load the raw data 
uniprot_data = pd.read_csv('../../data/schmidt2016_uniprot.tab', delimiter='\t')

#%%
# Rename some columns to be more meaningful.
uniprot_data.rename(columns={'Interacts with': 'interactions', 
                     'Gene ontology (biological process)':'go_process',
                     'Gene ontology (molecular function)':'go_function',
                     'Cofactor':'cofactor'}, inplace=True)

# Generate an empty dataframes for biological process and cofactors
cycle_df = pd.DataFrame([])

# Group by protein and split each entry as appropriate. 
for g, d in uniprot_data.groupby(['search']):
    # Get the list of processes 
    if g in gene_name_dict.keys():
        processes = d['go_process'].values[0]
        gene_name = gene_name_dict[g]
        if str(processes) == 'nan': 
           cycle_dict = {'uniprot':g,
                         'gene':gene_name,
                         'cog_class':cog_class_dict[g],
                         'cog_desc':cog_desc_dict[g],
                         'cog_letter':cog_letter_dict[g],
                         'desc':gene_desc_dict[g],
                         'go_number': np.nan,
                         'pathway': 'no associated pathway'
                         }
           cycle_df  = cycle_df.append(cycle_dict, ignore_index=True)

        else: 
            processes = processes.split(';')
            # Iterate through each entry and append to the appropriate dataframe
            for i, process in enumerate(processes):
                # Split by the GO ID no
                pathway, go_no = process.split('[GO:')
                go_no = int(go_no[:-1])

                # Check for an extraneous spaces
                if pathway[0] == ' ':
                    pathway = pathway[1:]
                if pathway[-1] == ' ':
                    pathway = pathway[:-1]

                # Append information to the dataframes
                cycle_dict = {'uniprot':g,
                              'gene':gene_name,
                              'cog_class':cog_class_dict[g],
                              'cog_desc':cog_desc_dict[g],
                              'cog_letter':cog_letter_dict[g],
                              'desc':gene_desc_dict[g],
                              'go_number':go_no,
                              'pathway':pathway,
                              }
                cycle_df = cycle_df.append(cycle_dict, ignore_index=True)
cycle_df.to_csv('../../data/uniprot_pathway_associations.csv', index=False)
# %%
