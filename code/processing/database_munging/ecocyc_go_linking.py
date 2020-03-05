#%%
import numpy as np
import pandas as pd
import glob
import tqdm

# Load the ecocyc data and restrict only to useful information
genes = pd.read_csv('../../../data/ecocyc_raw_data/ecocyc_ecoli_genes.tab', delimiter='\t')
genes.rename(columns={'Gene Name':'gene_name', 
                      'Accession-1':'b_number', 
                      'Product':'gene_product',
                      'GO terms (biological process)':'go_process',
                      'GO terms (cellular component)': 'go_component',
                      'GO terms (molecular function)' : 'go_function'},
                      inplace=True)
genes['gene_name'] = genes['gene_name'].str.lower()
genes.drop(labels=['Component-Of', 'Left-End-Position', 'Right-End-Position'], axis=1, inplace=True)

# Load the COG information
cogs = pd.read_csv('../../../data/ecoli_cog_gene_list.csv')
cogs['gene_name'] = cogs['gene_name'].str.lower()
cogs = cogs[['gene_name', 'cog_letter']]

#%%
# Manually flatten this into a dataframe where each row corresponds to a single
# go number. 
go_genes_df = pd.DataFrame([])
for g, d in tqdm.tqdm(genes.groupby(['gene_name', 'b_number', 'gene_product'])):
    annotation = 0
    for go in ['go_process', 'go_component', 'go_function']:
        processes = d[go].values[0]
        if str(processes) != 'nan':
            _processes = processes.split('//')
            for i, p in enumerate(_processes):
                # Grab the COG information, if it exists. 
                _cog = cogs[cogs['gene_name']==g[0]]
                if len(_cog) != 0:
                    letter = _cog['cog_letter'].values[0]
                    if len(letter) > 1:
                     letter = ''.join(list(letter))
                else:
                    letter = 'No Annotation'
                if i == 0:
                    _p = p.strip()
                else:
                    _p = p.strip()
                    _p = _p.replace('"', '')
                go_genes_df = go_genes_df.append({
                        'gene_name': g[0],
                        'b_number':g[1],
                        'gene_product':g[2],
                        'go_term': _p,
                        'go_class': go.split('_')[1],
                        'cog_letter': letter
                        }, ignore_index=True)   
        else:
            # Grab the COG information, if it exists. 
            _cog = cogs[cogs['gene_name']==g[0]]
            if len(_cog) != 0:
                letter = _cog['cog_letter'].values[0]
                if len(letter) > 1:
                    letter = ''.join(list(letter))
            else:
                letter = 'No Annotation'
            go_genes_df = go_genes_df.append({
                        'gene_name': g[0],
                        'b_number':g[1],
                        'gene_product':g[2],
                        'go_term': 'No Ontology',
                        'go_class': go.split('_')[1],
                        'cog_letter': letter
                        }, ignore_index=True)

#%% Link the cog letters to the appropriate desc and class
schmidt_data = pd.read_csv('../../../data/schmidt2016_longform.csv')
for g, d in schmidt_data.groupby(['cog_class', 'cog_desc', 'cog_class_letter']):
    go_genes_df.loc[go_genes_df['cog_letter']==g[-1], 'cog_class'] = g[0]
    go_genes_df.loc[go_genes_df['cog_letter']==g[-1], 'cog_category'] = g[1]

go_genes_df.loc[go_genes_df['cog_letter']=='No Annotation', 'cog_class'] = 'No Annotation'
go_genes_df.loc[go_genes_df['cog_letter']=='No Annotation', 'cog_category'] = 'No Annotation'
#%%
go_genes_df.to_csv('../../../data/ecoli_gene_list_go_cog.csv', index=False)


# %%
