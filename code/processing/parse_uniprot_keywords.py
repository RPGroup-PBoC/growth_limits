#%%
import numpy as np
import pandas as pd

# Define the identifiers
line_identifiers = {
    'ID': 'id',
    'AC': 'accession_number',
    'GO': 'go_mapping',
    'HI': 'hierarchy',
    'CA': 'category'
}

# Open the file as plaintext
with open('../../data/uniprot_kwlist.txt', 'r') as file:
    kwlist = file.read()

# Split the individual keywords
entries = kwlist.split('//\n')

df = pd.DataFrame([])
for entry in entries[:-1]:
    # Split the entries by lines:
    lines = entry.split('\n')

    # Find the important values
    kw_dict = {}
    for l in lines:
        for k, v in line_identifiers.items():
            if l[:2] == k:
                if k == 'HI':
                    kw_dict[v] = l.split(':')[1].split('.')[0]
                else:
                    kw_dict[v] = l.split(f'{k}   ')[1].split('.')[0]
    df = df.append(kw_dict, ignore_index=True)
df.to_csv('../../data/uniprot_keywords.csv', index=False)
# %%
