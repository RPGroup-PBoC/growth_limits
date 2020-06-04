#%%
import numpy as np
import pandas as pd
import tqdm
from scipy import constants
import prot.size as size

# Load the longform schmidt 2016 count data, growth rates, and cog associations
counts = pd.read_csv('../../../data/schmidt2016_raw_data/schmidt2016_dataset.csv')
rates = pd.read_csv('../../../data/schmidt2016_raw_data/schmidt2016_growth_rates.csv')
colicogs = pd.read_csv('../../../data/ecoli_genelist_master.csv')

# %%
# Get a list of the unique conditions.
conditions = rates['condition'].unique()

# Set up an empty dataframe
dfs = []

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby('gene'), desc="Iterating through genes..."):

    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name'].str.lower()==g.lower()]
    if len(gene) > 0:
        cog_class = gene['cog_class'].unique()[0]
        cog_cat = gene['cog_category'].unique()[0]
        cog_letter = gene['cog_letter'].unique()[0]
        b_number = gene['b_number'].unique()[0]
        gene_product = gene['gene_product'].unique()[0]
        go_term = ';'.join(list(gene['go_terms'].unique()))
        mw = gene['mw_fg'].values[0]

        # Iterate through each condition and extract relevant information.
        for c in conditions:
            growth_rate = rates.loc[rates['condition']==c]['growth_rate_hr'].values[0]
            reported_volume = rates.loc[rates['condition']==c]['volume_fL'].values[0]
            gene_dict = {
                    'gene_name': g,
                    'b_number': b_number,
                    'condition': c,
                    'reported_tot_per_cell': d[f'{c}_tot'].values[0],
                    'reported_fg_per_cell': d[f'{c}_tot'].values[0] * mw,
                    'go_terms': go_term,
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'growth_rate_hr': growth_rate,
                    'gene_product': gene_product,
                    'reported_volume': reported_volume,
                    'corrected_volume': size.lambda2size(growth_rate)
                    }
            dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(f'Warning!!! {g} not found in the gene list!')
        for c in conditions:
            growth_rate = rates.loc[rates['condition']==c]['growth_rate_hr'].values[0]
            reported_volume = rates.loc[rates['condition']==c]['volume_fL'].values[0]
            gene_dict = {
                    'gene_name': g[0],
                    'b_number': g[1],
                    'condition': c,
                    'reported_tot_per_cell': d[f'{c}_tot'].values[0],
                    'reported_fg_per_cell': d[f'{c}_fg'].values[0],
                    'cog_class': 'Not Found',
                    'cog_category': 'Not Found',
                    'cog_letter': 'Not Found',
                    'go_terms': 'Not Found',
                    'growth_rate_hr': growth_rate,
                    'reported_volume': reported_volume,
                    'corrected_volume': size.lambda2size(growth_rate)
                    }
            dfs.append(pd.DataFrame(gene_dict, index=[0]))
df = pd.concat(dfs, sort=False)

#%%
# Iterate through the conditions and correct fg and tot_per_cell as necessary.
for g, d in rates.groupby(['condition']):
    rel_corr_fg = df[df['condition'] == g]['reported_fg_per_cell'].sum() / \
                    size.lambda2P(d['growth_rate_hr'].unique())
                    
    print(g, ': total mass fg: ', np.round(size.lambda2P(d['growth_rate_hr'].unique()),2),
            ' volume: ', np.round(size.lambda2size(d['growth_rate_hr'].unique()),2),
            ' relative change in total fg: ', np.round(1/rel_corr_fg,2))

    df.loc[df['condition']==g, 'tot_per_cell'] = df.loc[df['condition']==g]['reported_tot_per_cell'] / rel_corr_fg
    df.loc[df['condition']==g, 'fg_per_cell'] =  df.loc[df['condition']==g]['reported_fg_per_cell'] / rel_corr_fg

#%%
df['dataset'] = 'schmidt_2016'
df['dataset_name'] = 'Schmidt et al. 2016'
df['strain'] = 'BW25113'
df.to_csv('../../../data/schmidt2016_longform_annotated.csv', index=False)
