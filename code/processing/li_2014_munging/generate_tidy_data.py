#%%
import numpy as np
import pandas as pd
import tqdm
# Define teh constant of total number of genes identified. 
N_GENES = 3041

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
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')

# Instantiate a blank dataframe to which measurements and annotations will be
# added
df = pd.DataFrame([])
for g, d in tqdm.tqdm(synthesis_tidy.groupby('Gene'), desc='Iterating through genes'):
    gene = colicogs.loc[colicogs['gene_name']==g]

    if len(g) > 0:
        for i in range(len(gene)):
            _gene = gene.iloc[i]
            cog_class = _gene['cog_class']
            cog_cat = _gene['cog_category']
            cog_letter = _gene['cog_letter']
            mass = _gene['mass_da']
            annotation = _gene['annotation']

            # Iterate through each condition and extract relevant information. 
            for c, r in growth_rates.items():
                if c in d['variable'].unique():
                    gene_dict = {
                        'gene_name': g,
                        'condition': c,
                        'tot_per_cell': d[d['variable']==c]['value'].values[0],
                        'cog_class': cog_class,
                        'cog_category': cog_cat,
                        'cog_letter': cog_letter,
                        'mass_da': mass,
                        'annotation': annotation,
                        'growth_rate_hr': growth_rates[c]
                        }
                    df = df.append(gene_dict, ignore_index=True)

#%%
# Compute the mass per cell and include dataset notation. 
df['fg_per_cell'] = df['mass_da'].values * df['tot_per_cell'].values * 6.022E-8
df['dataset'] = 'li_2014'
df['strain'] = 'MG1655'
df.to_csv('../../../data/li2014_longform_annotated.csv')
# %%
