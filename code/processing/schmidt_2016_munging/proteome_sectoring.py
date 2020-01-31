#%%
import pandas as pd
import prot.stats
from tqdm import tqdm

# Load the dataset (s)
condition_data = pd.read_csv('../../data/schmidt2016_longform.csv')
genes = pd.read_csv('../../data/schmidt2016_genes_processes.csv')

#%%
# Compute the total fraction occupied at each growth rate
cog_class_dfs = []
cog_desc_dfs = []
cog_gene_dfs = []
for g, d in tqdm(condition_data.groupby(['condition', 'growth_rate_hr']), 
    desc="Iterating through conditions"):
    # Find the global fractionation
    cog_class_frac = prot.stats.compute_fraction(d, 'cog_class')
    cog_class_frac['condition'] = g[0]
    cog_class_frac['growth_rate_hr'] = g[1]
    cog_class_dfs.append(cog_class_frac)

    # Find the subclass fractionation
    for _g, _d in tqdm(d.groupby('cog_class'), 
        desc="Iterating through cog classes"):
        cog_desc_frac = prot.stats.compute_fraction(_d, 'cog_desc')
        cog_desc_frac['condition'] = g[0]
        cog_desc_frac['growth_rate_hr'] = g[1]
        cog_desc_frac['cog_class'] = _g
        cog_desc_dfs.append(cog_desc_frac)

        # Find the gene fractionation
        for __g, __d in tqdm(_d.groupby('cog_desc'),
            desc="Iterating through genes"):
            cog_gene_frac = prot.stats.compute_fraction(__d, 'gene')
            cog_gene_frac['condition'] = g[0]
            cog_gene_frac['growth_rate_hr'] = g[1]
            cog_gene_frac['cog_class'] = _g
            cog_gene_frac['cog_desc'] = __g
            cog_gene_frac['desc'] = __d['desc']

            # ITerate through each gene, get the description, and prune. 
            for ___g, ___d in __d.groupby('gene'):
                _desc = ___d['desc'].unique()[0].split('OS')[0]
                cog_gene_frac.loc[cog_gene_frac['group']==___g, 'desc'] = _desc
            cog_gene_dfs.append(cog_gene_frac)

cog_class_df = pd.concat(cog_class_dfs, sort=False)
cog_class_df.to_csv('../../data/schmidt2016_cog_class_sectoring.csv', index=False)
cog_desc_df = pd.concat(cog_desc_dfs, sort=False)
cog_desc_df.to_csv('../../data/schmidt2016_cog_desc_sectoring.csv', index=False)
cog_gene_df = pd.concat(cog_gene_dfs, sort=False)
cog_gene_df.to_csv('../../data/schmidt2016_cog_gene_sectoring.csv', index=False)

#%%
# Compute the fractionation of each uniprot biological process by gene
process_gene_dfs = []
for g, d in tqdm(genes.groupby(['condition', 'uniprot_bio_process']),
            desc='Fractioning biological process gene constituents'):
    # Compute the fraction. 
    frac = prot.stats.compute_fraction(d, 'gene')
    frac['condition'] = g[0]
    frac['uniprot_bio_process'] = g[1]
    frac['growth_rate_hr'] = d['growth_rate_hr'].values[0]
    for _g, _d in frac.groupby('group'):
        frac.loc[frac['group']==_g, 'desc'] = d[d['gene']==_g]['desc'].values[0]
        frac.loc[frac['group']==_g, 'cog_class'] = d[d['gene']==_g]['cog_class'].values[0]
    process_gene_dfs.append(frac)

process_gene_df = pd.concat(process_gene_dfs)
process_gene_df.to_csv('../../data/schmidt2016_uniprot_process_gene_sectoring.csv', 
                        index=False)

#%%
# Compute fractionation of processes writ large
process_df = pd.DataFrame([])
for g, d in genes.groupby(['condition', 'uniprot_bio_process']):
    # Get the total mass and size of the proteome in this condition. 
    proteome_mass = condition_data[
                    condition_data['condition']==g[0]]['fg_per_cell'].sum()
    proteome_size = condition_data[
                    condition_data['condition']==g[0]]['tot_per_cell'].sum()

    # Compute the fractioning
    mass_frac = d['fg_per_cell'].sum() / proteome_mass
    count_frac = d['tot_per_cell'].sum() / proteome_size

    # Append information to the dataframe. 
    process_dict = {'condition':g[0], 
                    'uniprot_bio_process':g[1],
                    'growth_rate_hr': d['growth_rate_hr'].values[0],
                    'frac_mass': mass_frac,
                    'frac_count': count_frac}
    process_df = process_df.append(process_dict, ignore_index=True)

process_df.to_csv('../../data/schmidt2016_uniprot_process_sectoring.csv', index=False)

# %%
