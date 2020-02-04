#%%
import numpy as np
import pandas as pd
import tqdm

# Load the count data.
counts = pd.read_csv('../../../data/tanguichi2010_raw_data/tanguichi2010_counts.csv')
t_double = 150 # In min
# Load the cog annotation list.
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')

# Create the data frame using only the mean proteoin counts for now.
df = pd.DataFrame([])
for g, d in tqdm.tqdm(counts.groupby('Gene Name'), desc='Iterating through genes'):
    gene = colicogs.loc[colicogs['gene_name']==g]

    if len(g) > 0:
        for i in range(len(gene)):
            _gene = gene.iloc[i]
            cog_class = _gene['cog_class']
            cog_cat = _gene['cog_category']
            cog_letter = _gene['cog_letter']
            mass = _gene['mass_da']
            annotation = _gene['annotation']
            gene_dict = {
                        'gene_name': g,
                        'condition': 'M9_minimal_glucose_30C',
                        'tot_per_cell': d['Mean_Protein'].values[0],
                        'cog_class': cog_class,
                        'cog_category': cog_cat,
                        'cog_letter': cog_letter,
                        'mass_da': mass,
                        'annotation': annotation,
                        'growth_rate_hr': np.log(2) / (t_double / 60)
                        }
            df = df.append(gene_dict, ignore_index=True)

# Compute the mass per cell and define constants.
df['fg_per_cell'] = df['tot_per_cell'].values * df['mass_da'].values * 6.8E-8
df['dataset'] = 'tanguichi_2010'
df['strain'] = 'unknown'
df.to_csv('../../../data/tanguichi2010_longform_annotated.csv', index=False)

# %%
