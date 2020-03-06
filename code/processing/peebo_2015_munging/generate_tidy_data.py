#%%
import numpy as np
import pandas as pd
import tqdm

# Load the data quantifying absolute protein synthesis rates.
counts = pd.read_csv('../../../data/peebo2015_raw_data/peebo2014_copynums_minimal.csv')
counts = counts.append(pd.read_csv('../../../data/peebo2015_raw_data/peebo2014_copynums_rich.csv'))

# Get a list of the unique conditions.
conditions = counts['condition'].unique()[:2]

# Load the annotation list.
colicogs = pd.read_csv('../../../data/ecoli_genelist_master.csv')


# Instantiate a blank dataframe to which measurements and annotations will be
# added
dfs = []

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby('gene'), desc="Iterating through genes..."):

    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name'].str.lower()==g.lower()]
    b_number = gene['b_number'].unique()[0]
    gene_product = gene['gene_product'].unique()[0]
    go_term = ';'.join(list(gene['go_terms'].unique()))
    if len(gene) > 0:
        cog_class = gene['cog_class'].values[0]
        cog_cat = gene['cog_category'].values[0]
        cog_letter = gene['cog_letter'].values[0]
        gene_product= gene['gene_product'].values[0]
        mw = gene['mw_fg'].values[0]
        go_term = ';'.join(list(gene['go_terms'].unique()))
        for _c, _d in d.groupby(['growth_rate_hr-1']):
            vol = 0.27*2**(0.76*_c)
            # extract relevant information.
            gene_dict = {
                'gene_name': g.lower(),
                'b_number': b_number,
                'condition': _d['condition'].unique()[0],
                'corrected_volume': vol,
                'reported_tot_per_cell': _d['copy_number_molecule-per-fL'].values[0] * vol,
                'reported_fg_per_cell': _d['copy_number_molecule-per-fL'].values[0] * vol * mw,
                'go_terms':go_term,
                'cog_class': cog_class,
                'cog_category': cog_cat,
                'cog_letter': cog_letter,
                'gene_product': gene_product,
                'growth_rate_hr': _c
                }
            dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(f'Warning!!! {g} not found in the gene list! Not including in final tally')


#%%
# Compute the mass per cell and include dataset notation.
df = pd.concat(dfs, sort=False)
df['dataset'] = 'peebo_2015'
df['dataset_name'] = 'Peebo et al. 2015'
df['strain'] = 'BW25113'

#%%
# Compute the volume corrections 
_conditions = df.groupby(['growth_rate_hr', 'corrected_volume']).sum().reset_index()
_conditions['concentration'] = _conditions['reported_fg_per_cell'].values / _conditions['corrected_volume']
# Compute the relative concentration. 
rel_conc = _conditions[_conditions['growth_rate_hr']==0.55]['concentration'].values[0]
_conditions['rel_conc_to_ref'] = _conditions['concentration'] / rel_conc

#%% Update the counts. 
for g, d in _conditions.groupby(['growth_rate_hr']):
    rel_conc = _conditions[_conditions['growth_rate_hr']==g]['rel_conc_to_ref'].values[0]
    df.loc[df['growth_rate_hr']==g, 'tot_per_cell'] = df.loc[df['growth_rate_hr']==g]['reported_tot_per_cell'] / rel_conc
    df.loc[df['growth_rate_hr']==g, 'fg_per_cell'] =  df.loc[df['growth_rate_hr']==g]['reported_fg_per_cell'] / rel_conc


#%%
df.to_csv('../../../data/peebo2015_longform_annotated.csv')
# %%


# %%
