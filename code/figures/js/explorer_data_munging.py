#%%
import numpy as np
import pandas as pd
import tqdm
import prot.viz
dataset_colors = prot.viz.dataset_colors()
# ##############################################################################
# DATA CLEANING AND AGGREGATION
# ##############################################################################
# Load the three datasets
prots = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
cplx = pd.read_csv('../../../data/compiled_annotated_complexes.csv')
cplx = cplx[['gene_name', 'b_number', 'condition', 'growth_rate_hr', 'go_terms',
             'complex', 'complex_annotation', 'dataset', 'dataset_name',
             'n_subunits', 'n_units', 'gene_product']]
cplx = cplx[cplx['complex_annotation'] != 'none assigned']
cplx = cplx[cplx['complex'] != 'none assigned']
cplx

# Condense the complex data frame to an easily paresable form. 
cplx_numeric_dfs = []
cplx_desc_dfs = []

for g, d in tqdm.tqdm(cplx.groupby(['complex_annotation', 'complex']), 
                                    desc='Condensing complex data sets...'):

    cplx_desc = pd.DataFrame([])
    for _g, _d in d.groupby(['gene_name', 'n_subunits', 'gene_product']):
        formula += f' [{_g[0][0].upper()}{_g[0][1:]}]<sub>{int(_g[1])}</sub>'

        # Assemble a descriptive data frame
           cplx_desc = cplx_desc.append({'complex_annotation':g[0],
                                  'complex': g[1],
                                  'protein': _g[0][0].upper() + _g[0][1:],
                                  'subunits': _g[1], 
                                  'func': _g[2]},
                                  ignore_index=True)

    cplx_desc_dfs.append(cplx_desc)
    # Iterate through each complex, dataset, and condition, and compute the aggs
    cplx_numeric = pd.DataFrame([])
    for _g, _d in d.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
        cplx_numeric = cplx_numeric.append({
                          'min':_d['n_units'].min(),
                          'max':_d['n_units'].max(),
                          'mean':_d['n_units'].mean(),
                          'median':_d['n_units'].median(),
                          'complex': g[1],
                          'dataset': _g[0],
                          'dataset_name':_g[1],
                          'condition': _g[2],
                          'growth_rate_hr':_g[3],
                          'color': dataset_colors[_g[0]]
                          }, ignore_index=True)
    cplx_numeric_dfs.append(cplx_numeric)


# Define the column data sources
cplx_desc = pd.concat(cplx_desc_dfs, sort=False)
cplx_numeric = pd.concat(cplx_numeric_dfs, sort=False)
cplx_desc.to_csv('./cplx_desc.csv', index=False)
cplx_numeric.to_csv('./cplx_numeric.csv', index=False)

# %%
