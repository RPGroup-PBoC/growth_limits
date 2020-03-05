#%%
import numpy as np
import pandas as pd
import tqdm

# Load the data quantifying absolute protein synthesis rates.
counts = pd.read_csv('../../../data/valgepea2013_raw_data/valgepea2013_copynums.csv')

# Load the annotation list.
colicogs = pd.read_csv('../../../data/ecoli_genelist_master.csv')

# Instantiate a blank dataframe to which measurements and annotations will be
# added
dfs = []

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby(['gene', 'growth_rate_hr-1']), desc="Iterating through genes..."):
    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name'].str.lower()==g[0].split(',')[0].lower()]
      
    if len(gene) > 0:
        cog_class = _gene['cog_class']
        cog_cat = _gene['cog_category']
        cog_letter = _gene['cog_letter']
        gene_product = _gene['gene_product']
        go_term = ';'.join(list(gene['go_terms'].unique()))
        mw = gene['mw_kda'].values[0]
        if str(mw) == 'nan':
            mw = str(0)
        elif '//' in str(mw):
            mw = mw.split('//')[0]
        mw = float(mw)
        # extract relevant information.
        gene_dict = {
                    'gene_name': gene['gene_name'].unique()[0],
                    'b_number':b_number,
                    'condition': 'glucose_minimal',
                    'tot_per_cell': d['copy_number'].values[0],
                    'fg_per_cell':d['copy_number'].values[0] * mw * 1E3 * 6.022E-8, 
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'annotation': gene_product,
                    'growth_rate_hr': g[1],
                    'go_terms':go_term,
        }
        dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(g[0].split(',')[0])

#%%
# Compute the mass per cell and include dataset notation.
df['dataset'] = 'valgepea_2013'
df['dataset_name'] = 'Valgepea et al. 2013'
df['strain'] = 'MG1655'
df.to_csv('../../../data/valgepea2013_longform_annotated.csv')
# %%
