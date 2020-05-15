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
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
# %%
ecocyc = pythoncyc.select_organism('ECOLI')

prot_cplx = list(ecocyc.all_protein_complexes())
for frame in ecocyc['Protein-RNA-Complexes']:
    prot_cplx.append(frame['frameid'])

#%%
# Instantiate the dataframe.
_df = pd.DataFrame([])
for c in tqdm.tqdm(prot_cplx):
    components, counts = ecocyc.monomers_of_protein(c)
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

            # Get the annotation
            if type(ecocyc[c]['common_name']) != str:
                common_name = ecocyc[ecocyc[c]['catalyzes'][0]]['common_name']
            else:
                common_name = ecocyc[c]['common_name']

            # Update the dataframe
            _df = _df.append({'gene_name':gene.lower(), 
                              'n_copies':int(sub), 
                              'b_number': num,
                              'complex':c.strip('|'),
                              'annotation': common_name}, ignore_index=True)
    else:
        print(f'Barfed on {c}')

_dfs = []
for g, d in data.groupby(['gene_name', 'b_number', 'gene_product']):
    if g[1] not in list(_df['b_number'].values):
      _df = _df.append({'gene_name': g[0],
                         'n_copies': 1.0, 
                         'b_number': g[1],
                         'complex': 'none assigned',
                         'annotation': g[2]},
                         ignore_index=True)

_df.to_csv('../../../data/ecocyc_raw_data/annotated_complexes.csv', index=False)



# %%
