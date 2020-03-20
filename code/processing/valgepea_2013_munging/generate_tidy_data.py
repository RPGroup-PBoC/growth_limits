#%%
import numpy as np
import pandas as pd
import tqdm

# Load the data quantifying absolute protein synthesis rates.
counts = pd.read_csv('../../../data/valgepea2013_raw_data/valgepea2013_copynums.csv')
counts.dropna(axis=0, inplace=True)

#%%
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
        cog_class = gene['cog_class'].values[0]
        cog_cat = gene['cog_category'].values[0]
        cog_letter = gene['cog_letter'].values[0]
        gene_product = gene['gene_product'].values[0]
        go_term = ';'.join(list(gene['go_terms'].unique()))
        mw = gene['mw_fg'].values[0]
        b_number = gene['b_number'].values[0]
        # volume prediction Si, F. et al. (2017), Current Biology, http://doi.org/10.1016/j.cub.2017.03.022
        vol = 0.28 * np.exp(1.33  * g[1])
        # extract relevant information.
        gene_dict = {
                    'gene_name': gene['gene_name'].unique()[0],
                    'b_number':b_number,
                    'condition': 'glucose_minimal',
                    'reported_tot_per_cell': d['copy_number'].values[0],
                    'reported_fg_per_cell':d['copy_number'].values[0] * mw,
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'gene_product': gene_product,
                    'growth_rate_hr': g[1],
                    'go_terms':go_term,
                    'corrected_volume': vol
        }
        dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(f"Warning!!! Could not find {g[0].split(',')[0]} in gene list. Not including in final tally.")

# Compute the mass per cell and include dataset notation.
df = pd.concat(dfs, sort=False)
#%%
# Do the concentration correction
_conditions = df.groupby(['condition', 'growth_rate_hr', 'corrected_volume']).sum().reset_index()
_conditions['concentration'] = _conditions['reported_fg_per_cell'].values / _conditions['corrected_volume']
ref_conc = _conditions[_conditions['growth_rate_hr']==0.49]['concentration'].values[0]
_conditions['rel_conc_to_ref'] = _conditions['concentration'].values  / ref_conc

#%%
for g, d in _conditions.groupby(['growth_rate_hr']):
    rel_conc = _conditions[_conditions['growth_rate_hr']==g]['rel_conc_to_ref'].values[0]
    df.loc[df['growth_rate_hr']==g, 'tot_per_cell'] = df.loc[df['growth_rate_hr']==g]['reported_tot_per_cell'] / rel_conc
    df.loc[df['growth_rate_hr']==g, 'fg_per_cell'] =  df.loc[df['growth_rate_hr']==g]['reported_fg_per_cell'] / rel_conc

#%%
df['dataset'] = 'valgepea_2013'
df['dataset_name'] = 'Valgepea et al. 2013'
df['strain'] = 'MG1655'
df.to_csv('../../../data/valgepea2013_longform_annotated.csv')
# %%
