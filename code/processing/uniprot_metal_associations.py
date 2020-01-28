#%%
import pandas as pd 
import numpy as np

#Load the raw schmidt data to generate a dictionary of uniprot accession to name
data = pd.read_csv('../../data/schmidt2016_longform.csv')

# Load the raw data 
uniprot = pd.read_csv('../../data/schmidt2016_uniprot.tab', delimiter='\t')

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
uniprot.rename(columns={'Interacts with': 'interactions', 
                     'Gene ontology (biological process)':'go_process',
                     'Gene ontology (molecular function)':'go_function',
                     'Cofactor':'cofactor', 'Metal binding':'metal'}, inplace=True)

# Generate an empty dataframes for biological process and cofactors
metal_df = pd.DataFrame([])

# Group by protein and split each entry as appropriate. 
for g, d in uniprot.groupby(['searches']):
    # Deal with a mistake in the uniprot accession info
    if g == 'P00452,P00452-2':         
        g = 'P00452'
    if g == 'P02919,P02919-2':
        g = 'P02919'

    # Isolate the metal column. 
    metals = d['metal'].values[0]
    if str(metals) != 'nan':
        metal_coord_df = pd.DataFrame([])
        for entry in metals.split('METAL')[1:]:
            split = entry.split('/note="')[1].split(';')[0]
            metal_coord_df = metal_coord_df.append({'metal_id': split}, ignore_index=True)
    #         metal_dict = {
    #                      'metal':metal_id,
    #                      'coordination': count,
    #                      'uniprot':g,
    #                      'gene': gene_name_dict[g],
    #                      'desc': gene_desc_dict[g],
    #                      'cog_desc':  cog_desc_dict[g],
    #                      'cog_class': cog_class_dict[g],
    #                      'cog_class_letter': cog_letter_dict[g]}
    #         metal_df = metal_df.append(metal_dict, ignore_index=True)
    # else:
    #     metal_dict = {
    #                  'metal': 'none',
    #                  'coordination': 0,
    #                  'uniprot':g,
    #                  'gene': gene_name_dict[g],
    #                  'desc': gene_desc_dict[g],
    #                  'cog_desc':  cog_desc_dict[g],
    #                  'cog_class': cog_class_dict[g],
    #                  'cog_class_letter': cog_letter_dict[g]}
    #     metal_df = metal_df.append(metal_dict, ignore_index=True)



#%%

#     # Get the list of cofactors. 
#     cofactors = d['cofactor'].values[0]
#     if str(cofactors) != 'nan':
#         names = []
#         notes = []
#         for entry in cofactors.split('COFACTOR:')[1:]:
#             compounds = entry.split('Name=')[1:]
#             descriptions = entry.split('Note=')[1:]
#             cofactor = [c.split(';')[0] for c in compounds]
#             desc = [n.split(';')[0].split('{')[0][:-2] for n in descriptions]
#             for co, de in zip(cofactor, desc):
#                 # Assemble the dictionary
#                 cofactor_dict = {
#                                 'cofactor':co,
#                                 'cofactor_desc': de,
#                                 'uniprot':g,
#                                 'gene': gene_name_dict[g],
#                                 'desc': gene_desc_dict[g],
#                                 'cog_desc':  cog_desc_dict[g],
#                                 'cog_class': cog_class_dict[g],
#                                 'cog_class_letter': cog_letter_dict[g]
#                                 }
#                 cofactor_df = cofactor_df.append(cofactor_dict, ignore_index=True)
#     else:
#         cofactor_dict = {
#                         'cofactor':'None',
#                         'cofactor_desc': 'No known associated cofactor',
#                         'uniprot':g,
#                         'gene': gene_name_dict[g],
#                         'desc': gene_desc_dict[g],
#                         'cog_desc':  cog_desc_dict[g],
#                         'cog_class': cog_class_dict[g],
#                         'cog_class_letter': cog_letter_dict[g]
#                         }
#         cofactor_df = cofactor_df.append(cofactor_dict, ignore_index=True)
# cofactor_df.to_csv('../../data/uniprot_cofactors.csv', index=False)


# # %%


# %%
