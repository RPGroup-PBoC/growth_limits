#%%
import numpy as np
import pandas as pd
import tqdm

# Load the data quantifying absolute protein synthesis rates.
counts = pd.read_csv('../../../data/valgepea2013_raw_data/valgepea2013_copynums.csv')

# Load the annotation list.
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')

# Instantiate a blank dataframe to which measurements and annotations will be
# added
df = pd.DataFrame([])

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby(['gene', 'growth_rate_hr-1']), desc="Iterating through genes..."):

    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name'].str.lower()==g[0].split(',')[0]]
    if len(g) > 0:
        for i in range(len(gene)):
            _gene = gene.iloc[i]
            cog_class = _gene['cog_class']
            cog_cat = _gene['cog_category']
            cog_letter = _gene['cog_letter']
            mass = _gene['mass_da']
            annotation = _gene['annotation']
            # extract relevant information.
            gene_dict = {
                    'gene_name': gene['gene_name'].unique()[0],
                    'condition': 'glucose_minimal',
                    'tot_per_cell': d['copy_number'].values[0],
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'mass_da': mass,
                    'annotation': annotation,
                    'growth_rate_hr': g[1]
                    }
            df = df.append(gene_dict, ignore_index=True)

#%%
# Compute the mass per cell and include dataset notation.
df['fg_per_cell'] = df['mass_da'].values * df['tot_per_cell'].values * 6.022E-8
df['dataset'] = 'valgepea2013'
df['strain'] = 'MG1655'
df.to_csv('../../../data/valgepea2013_longform_annotated.csv')
# %%
