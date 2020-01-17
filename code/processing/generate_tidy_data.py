#%%
import pandas as pd
import tqdm

# Load the quantification measurements
quant = pd.read_csv('../../data/schmidt2016_dataset.csv')

# Load the COG assignments
cog = pd.read_csv('../../data/schmidt2016_cog_assignment.csv')

# Load the growth rate information
rates = pd.read_csv('../../data/schmidt2016_growth_rates.csv')

# %% 
# Generate sets of keys for renaming
conditions = [key.split('_fg')[0] for key in quant.keys() if '_fg' in key]

# Prune the descriptions
descs = [key.split('OS=')[0] for key in quant['desc'].values]

# Make an empty list for the data frames
dfs = []
for condition in conditions:
    _df = pd.DataFrame([])
    
    # Populate the dataframe with measured quantities
    _df['fg_per_cell'] = quant[condition + '_fg'].values
    _df['coeff_var'] = quant[condition + '_cv'].values
    _df['tot_per_cell'] = quant[condition + '_tot'].values

    # Populate the dataframe with identifying information
    _df['condition'] = condition
    _df['uniprot'] = quant['uniprot'].values
    _df['desc'] = descs
    _df['gene'] = quant['gene'].values
    _df['n_peptides'] = quant['n_peptides'].values
    _df['mw_da'] = quant['mw_da'].values

    # Find the growth rate for each condition.
    _df['growth_rate_hr'] = rates[
                    rates['condition']==condition]['growth_rate_hr'].values[0]
    
    # Append to the list and concatenate
    dfs.append(_df)
quant_longform = pd.concat(dfs)

# Iterate through each unique gene product and populate with the COG class information
for uniprot_id in quant_longform['uniprot'].unique():
    protein_cog = cog[cog['uniprot']==uniprot_id]
    quant_longform.loc[quant_longform['uniprot']==uniprot_id, 
                        'cog_class'] = protein_cog['cog_class'].values[0]
    quant_longform.loc[quant_longform['uniprot']==uniprot_id, 
                        'cog_desc'] = protein_cog['cog_desc'].values[0]
    quant_longform.loc[quant_longform['uniprot']==uniprot_id, 
                        'cog_class_letter'] = protein_cog['cog_desc'].values[0]

quant_longform.to_csv('../../data/schmidt2016_longform.csv', index=False)

# %%
