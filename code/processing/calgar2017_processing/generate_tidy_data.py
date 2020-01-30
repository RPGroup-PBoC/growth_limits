#%%
import numpy as np
import pandas as pd
import tqdm

# Load the metadata file and the complete protein count file. 
metadata = pd.read_csv('../../../data/caglar2017_raw_data/metaData.csv')
raw_data = pd.read_csv('../../../data/caglar2017_raw_data/proteinMatrix.csv')
raw_data = raw_data.melt('gene_id')
protein_id_keys = pd.read_csv('../../../data/caglar2017_raw_data/nameDictionary_RNA&Protein.csv')

# convert the protein_id_key to an easy lookup dictionary. 
protein_dict = {v[0][:-1]:v[1] for v in protein_id_keys[['Protein_id', 'gene_name']].values}

# Load the schmidt data for easy assignment of cog classes and categories. 
#%%
# Get the listed gene names
names = []
for g in raw_data['gene_id'].values:
    if g in protein_dict.keys():
        names.append(protein_dict[g])
    else:
        names.append('none')
raw_data['gene_name'] = names

#%%
# Drop all of the rows with nan in the gene name column. 
linked = raw_data[~raw_data['gene_name'].isnull()].copy()

# Replace the dataset name without the trailing value. 
dataset_names = [g[:-2] for g in linked['variable'].values]
linked['variable'] = dataset_names

#%%
# Iterate through each conditionand define the appropriate entries. 
iter = 0
for g in linked['variable'].unique():
    print(iter)
    iter += 1
    out =  metadata.loc[metadata['dataSet']==g][
                                ['experiment', 'growthTime_hr',
                                 'batchNumber', 'carbonSource', 
                                 'growthPhase', 'doublingTimeMinutes'
                                 ]]
    experiment, growth_time, batch_number, carbon_source, growth_phase, doubling_time = out.values[0]
    linked.loc[linked['variable']==g, 'experiment'] = experiment
    linked.loc[linked['variable']==g, 'carbon_source'] = carbon_source
    linked.loc[linked['variable']==g, 'growth_phase'] = growth_phase
    linked.loc[linked['variable']==g, 'growth_time_hr'] = growth_time
    linked.loc[linked['variable']==g, 'doubling_time'] = doubling_time

# Transform the doubling time to growth rate.
linked['growth_rate_hr'] = np.log(2) / linked['doubling_time'].values   
#%%
# Restrict the dataset solely to exponential phase growth. 
exponential = linked[linked['growth_phase']=='exponential']
#%%
exponential

# %%
