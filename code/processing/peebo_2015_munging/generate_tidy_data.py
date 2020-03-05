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
    if len(g) > 0:
        for i in range(len(gene)):
            _gene = gene.iloc[i]
            cog_class = _gene['cog_class']
            cog_cat = _gene['cog_category']
            cog_letter = _gene['cog_letter']
            gene_product= _gene['gene_product']
            go_term = ';'.join(list(gene['go_terms'].unique()))
            for _c, _d in d.groupby(['condition']):
                growth_rate = _d['growth_rate_hr-1'].unique()[0]
                # extract relevant information.
                gene_dict = {
                    'gene_name': gene['gene_name'].unique()[0],
                    'b_number': b_number,
                    'condition': _d['condition'].unique()[0],
                    'tot_per_cell': _d['copy_number_molecule-per-fL'].values[0]*0.27*2**(0.76*growth_rate),
                    'go_terms':go_term,
                    'cog_class': cog_class,
                    'cog_category': cog_cat,
                    'cog_letter': cog_letter,
                    'gene_product': gene_product,
                    'growth_rate_hr': growth_rate
                    }
            dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(f'Warning!!! {g} not found in the gene list!')
        for _g, _d in d.groupby('condition'):
            growth_rate = _d['growth_rate_hr-1'].unique()[0]
            gene_dict = {
                    'gene_name': g,
                    'b_number': b_number,
                    'condition': d,
                    'tot_per_cell': _d['copy_number_molecule-per-fL'].values[0]*0.27*2**(0.76*growth_rate),
                    'go_terms': 'Not Found',
                    'cog_class': 'Not Found',
                    'cog_category': 'Not Found',
                    'cog_letter': 'Not Found',
                    'go_terms': 'Not Found',
                    'growth_rate_hr': growth_rate}
            dfs.append(pd.DataFrame(gene_dict, index=[0]))

#%%
# Compute the mass per cell and include dataset notation.
df = pd.concat(dfs, sort=False)
df['dataset'] = 'peebo_2015'
df['dataset_name'] = 'Peebo et al. 2015'
df['strain'] = 'BW25113'
df.to_csv('../../../data/peebo2015_longform_annotated.csv')
# %%
