#%%
import numpy as np
import pandas as pd
import tqdm

# Load the cog associations. 
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')

# Load the total protein counts. 
counts = pd.read_csv('../../../data/soufi2015_raw_data/soufi2015_tp3_counts.csv')

# Define some hardcoded values -- note that the growth rate is assumed standard
# for M9 glucose (0.5%) 
t_double = 60 # in min
growth_rate = np.log(2)  / (t_double / 60)

#%%
# Iterate through each gene in the count data. 
df = pd.DataFrame([])
for g, d in counts.groupby('Protein IDs'):

    # Determine number of entries per gene
    if len(g.split(';')) > 1:
        gene = []
        for pid in g.split(';'):
            if len(gene) == 0:
                _gene = colicogs[colicogs['uniprot']==pid]
                if len(_gene) > 0:
                    gene = _gene
    else:
        gene = colicogs[colicogs['uniprot']==g]

    if len(g) > 0:
        for i in range(len(gene)):
            _gene = gene.iloc[i]
            gene_name = _gene['gene_name']
            cog_class = _gene['cog_class']
            cog_cat = _gene['cog_category']
            cog_letter = _gene['cog_letter']
            mass = _gene['mass_da']
            annotation = _gene['annotation']

            #
    
            # Iterate through each condition and extract relevant information. 
            gene_dict = {
                        'gene_name': gene_name,
                        'condition': 'M9_glucose',
                        'tot_per_cell': d['CN (TP3) Rep2TP3'].values[0],
                        'cog_class': cog_class,
                        'cog_category': cog_cat,
                        'cog_letter': cog_letter,
                        'mass_da': mass,
                        'annotation': annotation,
                        'growth_rate_hr': growth_rate
                        }
            df = df.append(gene_dict, ignore_index=True)

# Add iddnetifying information. 
df['fg_per_cell'] = df['tot_per_cell'].values * df['mass_da'].values * 6.8E-8
df['dataset'] = 'soufi_2015'
df['strain'] = 'BW25113'
df.to_csv('../../../data/soufi2015_longform_annotated.csv', index=False)
# %%
