#%%
import numpy  as np
import pandas as pd 
import tqdm
import glob
#%%
# Load the raw ecocyc master list. 
ecocyc = pd.read_csv('../../../data/ecocyc_raw_data/2020-03-04_ecocyc_master.tab',
                    delimiter='\t')

# Iterate through each gene name and, for each synonym, create a new entry. 
df = pd.DataFrame([])
for g, d in tqdm.tqdm(ecocyc.groupby(['gene_name', 'ecocyc_gene_id', 
                            'b_number', 'gene_product', 
                            'common_name', 'mw_kda'])):
    # Stitch together the go ids. 
    go_ids = d['go_terms'].values[0]
    if str(go_ids) != 'nan':
        go_ids = '; '.join([s.strip() for s in go_ids.split('//')])
    else:
        go_ids = 'no ontology'

    # Generate the base dict. 
    base_dict = {'gene_name':g[0].lower(),
                'ecocyc_gene_id':g[1],
                'b_number':g[2],
                'gene_product':g[3],
                'mw_kda':g[5],
                'go_terms':go_ids}

    # Update the dataframe
    df = df.append(base_dict, ignore_index=True)

    # Iterate through each synonym
    syn = d['synonyms'].values[0]
    if str(syn) != 'nan':
        syn_split = syn.split('//')
        syns = [s.strip().replace('"', '').lower() for s in syn_split]
        for s in syns:
            base_dict['gene_name'] = s
            df = df.append(base_dict, ignore_index=True)
#%%
# Drop duplicate rows
df.drop_duplicates(inplace=True)


# %%
# Load all of the cog lists and collate
cogs = pd.concat([pd.read_csv(f) for f in glob.glob('../../../data/cog_data/*.csv')])
cogs.drop(columns=['Unnamed: 4'], axis=1, inplace=True)

# Assign the cog information to the ecocyc df based on the b number
for g, d in df.groupby('b_number'):
    cog = cogs[cogs['b_number']==g]
    if len(cog) == 0:
        df.loc[df['b_number']==g, 'cog_class'] = 'Not Assigned'
        df.loc[df['b_number']==g, 'cog_letter'] = 'Not Assigned'
        df.loc[df['b_number']==g, 'cog_category'] = 'Not Assigned'
        df.loc[df['b_number']==g, 'cog_desc'] = 'Not Assigned'
    else:
        df.loc[df['b_number']==g, 'cog_class'] = cog['cog_class'].values[0]
        df.loc[df['b_number']==g, 'cog_letter'] = cog['cog_letter'].values[0]
        df.loc[df['b_number']==g, 'cog_category'] = cog['cog_category'].values[0]
        df.loc[df['b_number']==g, 'cog_desc'] = cog['cog_desc'].values[0]

# %%
df.to_csv('../../../data/ecoli_genelist_master.csv', index=False)

# %%
