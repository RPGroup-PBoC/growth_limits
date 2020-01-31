"""
This script associates all genes from *E. coli* K-12 to a variety of
identifiers, including RefSeq ID's, UniProt accession numbers, Protein Ids, and
COG subgroups adn COG classes.

"""
#%%
import pandas as pd

# Load the association of cog id to e coli gene. 
colicogs = pd.read_csv('../../../data/ecoli_cog_gene_list.csv')

# Load the dataset with the discontinued refseq gene associations
refseq_gene = pd.read_csv('../../../data/caglar2017_raw_data/protein_mRNA_dictionary.csv')
refseq_gene = refseq_gene[~refseq_gene['gene_name'].isnull()]

# Load the ecocyc database with ID mapping
ecocyc = pd.read_table('../../../data/ecocyc_20200101.tab.gaf-1', index_col=False)

# Prune the ecocyc list to relevant columns. 
ecocyc = ecocyc[['uniprot', 'name']]
ecocyc.rename(columns={'name':'gene_name'}, inplace=True)


# %% Munge some details of the numbering scheme
for g, d in ecocyc.groupby(['gene_name', 'uniprot']):
    colicogs.loc[colicogs['gene_name']==g[0], 'uniprot'] = g[1]

# %%
