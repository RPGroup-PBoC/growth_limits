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
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')

# cell volume esimate
# Taheri-Araghi et al. 2015 best fit: Y = 0.27*2$^{0.76*\lambda}
# vol = 0.27*2**(0.76*\lambda)

# Instantiate a blank dataframe to which measurements and annotations will be
# added
df = pd.DataFrame([])

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby(['gene', 'growth_rate_hr-1', 'condition']), desc="Iterating through genes..."):

    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name'].str.lower()==g[0]]
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
                    'condition': g[2],
                    'tot_per_cell': d[d['condition']==g[2]]['copy_number_molecule-per-fL'].values[0]*0.27*2**(0.76*g[1]),
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
df['dataset'] = 'peebo_2015'
df['strain'] = 'BW25113'
df.to_csv('../../../data/peebo2015_longform_annotated.csv')
# %%
