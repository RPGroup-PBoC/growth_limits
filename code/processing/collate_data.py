#%%
import pandas as pd
import glob

# Grab all of the longform data, read, concatenate, and save. 
files = glob.glob('../../../data/*longform_annotated.csv')
dfs = [pd.read_csv(f) for f in files]
df = pd.concat(dfs, sort=False)
df.to_csv('../../../data/compiled_datasets.csv', index=False)