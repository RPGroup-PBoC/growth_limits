"""
This script associates all genes from *E. coli* K-12 to a variety of
identifiers, including RefSeq ID's, UniProt accession numbers, Protein Ids, and
COG subgroups adn COG classes.

"""
#%%
import pandas as pd

# Load the tidy e coli protein list from Caglar2017
coli_prots = pd.read_csv('../../data/caglar2017_raw_data/protein_tidy_eColi_ez.csv',
            index_col='Unnamed: 0')

# Load the dataset with the discontinued refseq gene associations
refseq_gene = pd.read_csv('../../data/caglar2017_raw_data/protein_mRNA_dictionary.csv')

# Load the official (note, discontinued) RefSeq-ProteinID 
refseq_protid = pd.read_table('../../data/caglar2017_raw_data/RefSeq_to_ProtID.tab')

# Join the dataframes on the gene name
# %% Munge some details of the numbering scheme
