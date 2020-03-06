#%%
import numpy as np
import pandas as pd
import tqdm

# From the SI of Li 2014, Cell, hard-code the growth rates. 
growth_rates = {'MOPS complete': np.log(2) / (21.5 / 60), 
        'MOPS complete without methionine': np.log(2) / (26.5 / 60),
        'MOPS minimal': np.log(2) / (56.3 / 60)}

# Load the data quantifying absolute protein synthesis rates. 
synthesis_rates = pd.read_csv('../../../data/li2014_raw_data/li2014_synthesis_rates.csv')

#%%
# Reform the synthesis rates to tidy-format 
synthesis_tidy = synthesis_rates.melt('Gene')
synthesis_tidy = synthesis_tidy[~synthesis_tidy['value'].str.contains('\[')]

# Load the annotation list. 
colicogs = pd.read_csv('../../../data/ecoli_genelist_master.csv')

# Instantiate a blank dataframe to which measurements and annotations will be
# added
dfs = [] 
for g, d in tqdm.tqdm(synthesis_tidy.groupby('Gene'), desc='Iterating through genes'):
    if '+' in g:
        _gene = g.split('+')
        split_factor = 2
    else:
        _gene = [g]
        split_factor=1
    for _g in _gene:
        gene = colicogs.loc[colicogs['gene_name'].str.lower()==_g.lower()]
        if len(gene) > 0:
            cog_class = gene['cog_class'].values[0]
            cog_cat = gene['cog_category'].values[0]
            cog_letter = gene['cog_letter'].values[0]
            b_number = gene['b_number'].unique()[0]
            gene_product = gene['gene_product'].unique()[0]
            go_term = ';'.join(list(gene['go_terms'].unique()))
            mw = gene['mw_fg'].values[0]
            # Iterate through each condition and extract relevant information. 
            for c, r in growth_rates.items():
                if c in d['variable'].unique():
                    gene_dict = {
                        'gene_name': _g,
                        'condition': c,
                        'tot_per_cell': split_factor**-1 * float(d[d['variable']==c]['value'].values[0]),
                        'fg_per_cell': split_factor**-1 * float(d[d['variable']==c]['value'].values[0]) * mw,
                        'cog_class': cog_class,
                        'cog_category': cog_cat,
                        'cog_letter': cog_letter,
                        'b_number':b_number,    
                        'gene_product': gene_product,
                        'growth_rate_hr': growth_rates[c],
                        'go_terms': go_term}
                    dfs.append(pd.DataFrame(gene_dict, index=[0]))
        else:
            print(f'Warning!! {g} not found in list.')
df = pd.concat(dfs, sort=False)
#%%
# Compute the mass per cell and include dataset notation. 
df['dataset'] = 'li_2014'
df['dataset_name'] = 'Li et al. 2014'
df['strain'] = 'MG1655'
df.to_csv('../../../data/li2014_longform_annotated.csv')
# %%
