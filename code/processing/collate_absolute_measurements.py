# %%
import pandas as pd
import glob

# Define the author names which did absolute measurments that are trustworth.
authors = ['schmidt2016', 'li2014', 'valgepea2013', 'tanguichi2010', 'peebo2015']

# Load inthe datasets and concatenate.
data = pd.concat([
    pd.read_csv(f'../../data/{auth}_longform_annotated.csv') for auth in authors],
    sort=False)
data.to_csv('../../data/compiled_absolute_measurements.csv', index=False)
