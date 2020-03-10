"""
Note that this script can only be run if BioCyc's pathway tools is installed 
locally and is being served via the following command:

./pathway-tools -lisp -python
"""
#%%
import numpy as np
import pandas as pd
import tqdm as tqdm
import pythoncyc

annotations = pd.read_csv('../../../data/ecoli_genelist_master.csv')
# %%
ecocyc = pythoncyc.select_organism('ECOLI')

# Get the list of protein complexes. 
cplx = ecocyc.all_protein_complexes()

# Instantiate the dataframe.
_df = pd.DataFrame([])
for c in tqdm.tqdm(cplx):
    components, counts = ecocyc.base_components_of_protein(c)
    if (components is not None) & (counts is not None):
        for prot, sub in zip(components, counts):
            if ecocyc[prot].abbrev_name is None:
                if ecocyc[prot].synonyms is None:
                    gene = prot.split('-')[0][1:]
                else:
                    gene = ecocyc[prot].synonyms
                    if type(gene) is list:
                        gene = gene[0]
            else:
                gene = ecocyc[prot].abbrev_name
                if type(gene) is list:
                    gene = gene[0]
            # Remove greek if necessary
            gene = gene.replace('&beta;', 'B')
            gene = gene.replace('&alpha;', 'A')

            # Get the b number
            num = annotations[annotations['gene_name'].str.lower() == gene.lower()]['b_number'].values

            if len(num) > 0:
                num = num[0] 

            # Update the dataframe
            _df = _df.append({'gene_name':gene.lower(), 
                              'n_copies':int(sub), 
                              'b_number': num,
                              'complex':c.strip('|'),
                              'annotation': ecocyc[c]['common_name']}, ignore_index=True)
    else:
        print(f'Barfed on {c}')

# %%
_df.to_csv('../../../data/ecocyc_raw_data/annotated_complexes.csv', index=False)



# %%
