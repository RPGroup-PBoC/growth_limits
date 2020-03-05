#%%
import numpy as np
import pandas as pd
import tqdm
from scipy import constants

# Load the longform schmidt 2016 count data, growth rates, and cog associations
counts = pd.read_csv('../../../data/schmidt2016_raw_data/schmidt2016_dataset.csv')
rates = pd.read_csv('../../../data/schmidt2016_raw_data/schmidt2016_growth_rates.csv')
colicogs = pd.read_csv('../../../data/ecoli_gene_list_go_cog.csv')

# %%
# Get a list of the unique conditions.
conditions = rates['condition'].unique()

# Set up an empty dataframe
df = pd.DataFrame([])

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby('gene'), desc="Iterating through genes..."):

    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name']==g.lower()]
    if len(gene) > 0:
        cog_class = gene['cog_class'].unique()[0]
        cog_cat = gene['cog_category'].unique()[0]
        cog_letter = gene['cog_letter'].unique()[0]
        b_number = gene['b_number'].unique()[0]
        gene_product = gene['gene_product'].unique()[0]
        go_term = ';'.join(list(gene['go_term'].unique()))
        go_class = ';'.join(list(gene['go_class'].unique()))
        # Iterate through each condition and extract relevant information.
        for c in conditions:
            gene_dict = {
                    'gene_name': g,
                    'condition': c,
                    'tot_per_cell': d[f'{c}_tot'].values[0],
                    'fg_per_cell':d[f'{c}_fg'].values[0],
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'growth_rate_hr': rates.loc[rates['condition']==c]['growth_rate_hr'].values[0]
                    }
            df = df.append(gene_dict, ignore_index=True)

# Compute the mass per cell.
# df['fg_per_cell'] = (df['tot_per_cell'].values * df['mass_da'] * 1E15 ) / constants.Avogadro 
df['dataset'] = 'schmidt_2016'
# df['strain'] = 'BW25153'
# df.to_csv('../../../data/schmidt2016_longform_annotated.csv', index=False)

# %%
