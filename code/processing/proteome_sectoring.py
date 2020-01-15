import pandas as pd
import prot.stats
from tqdm import tqdm

# Load the dataset 
data = pd.read_csv('../../data/schmidt2016_longform.csv')

# Compute the total fraction occupied at each growth rate
cog_class_dfs = []
cog_desc_dfs = []
cog_gene_dfs = []
for g, d in tqdm(data.groupby(['condition', 'growth_rate_hr']), 
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
            cog_gene_dfs.append(cog_gene_frac)

cog_class_df = pd.concat(cog_class_dfs, sort=False)
cog_class_df.to_csv('../../data/schmidt2016_cog_class_sectoring.csv', index=False)
cog_desc_df = pd.concat(cog_desc_dfs, sort=False)
cog_desc_df.to_csv('../../data/schmidt2016_cog_desc_sectoring.csv', index=False)
cog_gene_df = pd.concat(cog_gene_dfs, sort=False)
cog_gene_df.to_csv('../../data/schmidt2016_cog_gene_sectoring.csv', index=False)

