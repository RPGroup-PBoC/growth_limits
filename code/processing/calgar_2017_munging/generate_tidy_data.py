#%%
import numpy as np
import pandas as pd
import tqdm

# Load the metadata file and the complete protein count file. 
metadata = pd.read_csv('../../../data/caglar2017_raw_data/metaData.csv')
raw_data = pd.read_csv('../../../data/caglar2017_raw_data/proteinMatrix.csv')
raw_data = raw_data.melt('gene_id')
raw_data.rename(columns={'gene_id':'refseq_id', 
                         'variable':'experiment', 
                         'value':'tot_per_cell'},
                         inplace=True)

# Load the gene association file. 
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_annotations.csv')

# make a lookup table of refseq ids and gene name 
protein_names = {v[0]:v[1] for v in colicogs[['refseq_id', 'gene_name']].values}

# Convert the refseq id to the appropriate protein name and drop unindentified
# genes
names = []
for k in raw_data['refseq_id'].values:
    try:
        names.append(protein_names[k])
    except KeyError:
        names.append(np.nan)
raw_data['gene_name'] = names

# %% Prune the data to keep only genes present in the uniprot database
linked = raw_data[~raw_data['gene_name'].isnull()]

#%%
# Replace the dataset name without the trailing value. 
dataset_names = [g[:-2] for g in linked['experiment'].values]
linked['experiment'] = dataset_names.copy()

#%%
# Iterate through each condition and define the appropriate entries. 
for g in tqdm.tqdm(linked['experiment'].unique(), desc='Iterating through experiments'):
    # Isolate the corrected metadata entry
    out =  metadata.loc[metadata['dataSet']==g][
                                ['experiment', 'growthTime_hr',
                                 'batchNumber', 'carbonSource', 
                                 'growthPhase', 'doublingTimeMinutes'
                                 ]]
    experiment, growth_time, batch_number, carbon_source, growth_phase, doubling_time = out.values[0]
    
    # Add metadata information to the list
    linked.loc[linked['experiment']==g, 'condition'] = experiment
    linked.loc[linked['experiment']==g, 'carbon_source'] = carbon_source
    linked.loc[linked['experiment']==g, 'growth_phase'] = growth_phase
    linked.loc[linked['experiment']==g, 'growth_time_hr'] = growth_time
    linked.loc[linked['experiment']==g, 'doubling_time'] = doubling_time

# Transform the doubling time to growth rate.
linked['growth_rate_hr'] = np.log(2) / (linked['doubling_time'].values / 60)

#%%
# Assign the cog information to each gene and compute the mass per cell. 
dfs = []
annotated = linked.copy()
iter = 0
for g in annotated['gene_name'].unique():
    if iter %10 == 0:
        print(iter)
    iter += 1
    cog_info = colicogs.loc[colicogs['gene_name']==g]

    for i in range(len(cog_info)):
        _cog_info = cog_info.iloc[i].to_dict()
        if i == 0:
            annotated.loc[annotated['gene_name']==g, 
                                     'annotation'] = _cog_info['annotation']
            annotated.loc[annotated['gene_name']==g, 
                                     'cog_category'] = _cog_info['cog_category']
            annotated.loc[annotated['gene_name']==g, 
                                     'cog_class'] = _cog_info['cog_class']
            annotated.loc[annotated['gene_name']==g, 
                                     'cog_letter'] = _cog_info['cog_letter']
            annotated.loc[annotated['gene_name']==g, 
                                     'mass_da'] = _cog_info['mass_da']

        else:
            _d = annotated.loc[annotated['gene_name']==g].copy() 
            _d['annotation'] = _cog_info['annotation']
            _d['cog_category'] = _cog_info['cog_category']  
            _d['cog_class'] = _cog_info['cog_class']
            _d['cog_letter'] = _cog_info['cog_letter']
            _d['mass_da'] = _cog_info['mass_da']
            dfs.append(_d)

#%%    
# Concatenate the duplicate dataframes and fold into the annotation.
dup_copy = [d for d in dfs]
dup_copy.append(annotated)
complete_annotation = pd.concat(dup_copy, sort=False)

# Compute the mass per cell. 
complete_annotation['fg_per_cell'] = complete_annotation['tot_per_cell'].values *\
                                     complete_annotation['mass_da'].values * 6.022E-8

# Keep only the exponential phase growth at t = 6 hr. 
complete_annotation = complete_annotation[
    (complete_annotation['growth_phase']=='exponential') & 
    (complete_annotation['growth_time_hr']==6)]

# Save to disk. 
complete_annotation.drop(labels=['growth_phase', 'refseq_id'], axis=1, inplace=True)
complete_annotation['dataset'] = 'calgar_2017'
complete_annotation['strain'] = 'REL606'
complete_annotation.to_csv('../../../data/caglar2017_longform_annotated.csv', 
                    index=False)

# %%
