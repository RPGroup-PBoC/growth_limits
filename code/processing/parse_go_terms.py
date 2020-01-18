#%%
import pandas as pd
import tqdm

with open('../../data/go.obo', 'r') as file:
    lines = file.read()

# Split the gene ontology list by terms
terms = lines.split('[Term]')[1:]

# Set up the empty dataframe
df = pd.DataFrame([])
for term in tqdm.tqdm(terms, desc='Parsing GO terms...'):
    # Split by lines
    lines = term.split('\n')[1:]
    # Iterate through each line. 
    term_dict = {}
    for l in lines:
        if 'id:' == l[:3]:
            term_dict['id'] = l.split('id: ')[1]
        elif 'name:' == l[:5]:
            term_dict['name'] = l.split('name: ')[1]
        elif 'namespace:' == l[:10]:
            term_dict['namespace'] = l.split('namespace: ')[1]
        elif 'def:' == l[:5]:
            term_dict['definition'] = l.split('def: ')[1].split('[')[0][1:-1]
        if 'is_obsolete:' in l:
            term_dict['obsolete'] = True
        else: 
            term_dict['obsolete'] = False
        
    # Set up the new dictionary
    df = df.append(term_dict, ignore_index=True)                        
# %%
df.to_csv('../../data/go_terms.csv', index=False)


# %%
