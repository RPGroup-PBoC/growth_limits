#%%
import numpy as np
import pandas as pd
import tqdm

# From the SI of Li 2014, Cell, hard-code the growth rates.
growth_rate = np.log(2) / (150 / 60)

# Load the data quantifying absolute protein synthesis rates.
counts = pd.read_csv('../../../data/taniguchi2010_raw_data/taniguchi2010_counts.csv')
counts = counts[['gene_name','Mean_Protein', 'b_number']]
print(counts.head())
#%%
# Load the annotation list.
colicogs = pd.read_csv('../../../data/ecoli_genelist_master.csv')

# Instantiate a blank dataframe to which measurements and annotations will be
# added
dfs = []

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby(['b_number']), desc="Iterating through genes..."):
    # Determine number of entries per gene
    gene = colicogs[colicogs['b_number'].str.lower()==g]

    if len(gene) > 0:
        cog_class = gene['cog_class'].values[0]
        cog_cat = gene['cog_category'].values[0]
        cog_letter = gene['cog_letter'].values[0]
        gene_product = gene['gene_product'].values[0]
        go_term = ';'.join(list(gene['go_terms'].unique()))
        mw = gene['mw_fg'].values[0]
        b_number = g #gene['b_number'].values[0]
        # volume prediction Si, F. et al. (2017), Current Biology, http://doi.org/10.1016/j.cub.2017.03.022
        vol = 0.28 * np.exp(1.33  * growth_rate)
        # extract relevant information.
        gene_dict = {
                    'gene_name': gene['gene_name'].unique()[0],
                    'b_number':b_number,
                    'condition': 'M9_minimal_glucose_30C',
                    'reported_tot_per_cell': d['Mean_Protein'].values[0],
                    'reported_fg_per_cell':d['Mean_Protein'].values[0] * mw,
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'gene_product': gene_product,
                    'growth_rate_hr': growth_rate,
                    'go_terms':go_term,
                    'corrected_volume': vol
        }
        dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(f"Warning!!! Could not find {g} in gene list. Not including in final tally.")

# Compute the mass per cell and include dataset notation.
df = pd.concat(dfs, sort=False)

#%%
df['dataset'] = 'taniguichi_2010'
df['dataset_name'] = 'Taniguichi et al. 2010'
df['strain'] = 'DY330'
df.to_csv('../../../data/taniguchi2010_longform_annotated.csv')
# %%
