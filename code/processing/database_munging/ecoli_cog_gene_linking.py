#%%
import pandas as pd 
import tqdm 

#%% Load the schmidt data to set up a dictionary of cog letter to sub and super
#class
schmidt_data = pd.read_csv('../../../data/schmidt2016_longform.csv')
cog_info = {}
for g, d in schmidt_data.groupby(['cog_class', 'cog_desc', 'cog_class_letter']):
    cog_info[g[-1]] = g[1]
    cog_info[g[1]] = g[0]

#%%
# Load the file. 
with open('../../../data/COG_ID_org_prot.txt', 'r') as file:
    text = file.read()

# Split the parsed text by the entry distinctions
entries = text.split('_______')

#%%
# Set up an empty dataframe to hold the associations. 
colicogs = pd.DataFrame([])

# Iterate through each entry. 
for entry in tqdm.tqdm(entries[1:-1], desc="Processing COG entries"):
    subentries = entry.split('\n ')

    # Parse relevant ID information from the first line of the subentries. 
    COG_number = f"COG{subentries[0].split('COG')[1].split(' ')[0]}"
    annotation = subentries[0].split(f'{COG_number} ')[1]
    COG_letters = subentries[0].split('\n\n')[1].split('[')[1].split(']')[0]

    # Find the entry for ecoli (Eco)
    ecoli_entries = [s.split('Eco: ')[1] for s in subentries if 'Eco:' in s]    
    if len(ecoli_entries) > 0:
        for cog in COG_letters:
            for ec_entry in ecoli_entries:
                if len(ec_entry.split(' ')) > 1:
                    for gene in ec_entry.split(' '):
                        if gene != '': 
                            entry_dict = {'cog_id': COG_number,
                                  'cog_letter': cog,
                                  'cog_class': cog_info[cog],
                                  'cog_category': cog_info[cog_info[cog]],
                                  'gene_name': gene,
                                  'annotation': annotation}
                            colicogs = colicogs.append(entry_dict, ignore_index=True) 
                else:
                    entry_dict = {'cog_id': COG_number,
                                  'cog_letter': cog,
                                  'cog_class': cog_info[cog],
                                  'cog_category': cog_info[cog_info[cog]],
                                  'gene_name': ec_entry,
                                  'annotation': annotation}
                    colicogs = colicogs.append(entry_dict, ignore_index=True) 
colicogs = colicogs[~colicogs['gene_name'].isnull()]
#%%
# Load the dataset with the discontinued refseq gene associations
refseq_gene = pd.read_csv('../../../data/caglar2017_raw_data/protein_mRNA_dictionary.csv')
refseq_gene = refseq_gene[~refseq_gene['gene_name'].isnull()]

# Load the ecocyc database with ID mapping
ecocyc = pd.read_table('../../../data/ecocyc_20200101.tab.gaf-1', index_col=False)

#%%
ecocyc = ecocyc[['uniprot', 'name', 'go_id']]
ecocyc.rename(columns={'name':'gene_name'}, inplace=True)

# Load the uniprot entries with associated keywords for all ecoli genes. 
uniprot_kws = pd.read_table('../../../data/ecoli_uniprot_kws.tab')

# %% Assign uniprot and refseq id 
for g, d in tqdm.tqdm(ecocyc.groupby(['gene_name', 
                              'uniprot']), desc='Assigning UniProt Identifiers'):
    # Get the list of go_ids
    goids = '; '.join(list(d['go_id'].unique()))
    colicogs.loc[colicogs['gene_name']==g[0], 'uniprot'] = g[1]
    colicogs.loc[colicogs['gene_name']==g[0], 'go_ids'] = goids
    
for g, d in tqdm.tqdm(refseq_gene.groupby(['gene_name', 
                            'Protein_id']), desc='Assigning RefSeq Identifiers'):
    colicogs.loc[colicogs['gene_name']==g[0], 'refseq_id'] = g[1]

for g, d in tqdm.tqdm(uniprot_kws.groupby(['Entry', 'Keyword ID', 'Mass']),
                        desc='Assigning UniProt Keywords'):
    colicogs.loc[colicogs['uniprot']==g[0], 'uniprot_kw_ids'] = g[1]
    colicogs.loc[colicogs['uniprot']==g[0], 'mass_da'] = float(''.join(g[2].split(',')))


# Prune the dataframe to drop entries that are *not* present in the uniprot
# database (mostly an issue for Caglar2017 as of now 01/30/20)
colicogs = colicogs[~colicogs['uniprot'].isnull()]

# %% Save the final likage dataframe
colicogs.to_csv('../../../data/escherichia_coli_gene_annotations.csv', index=False)

# %%
