# %%
import pandas as pd
import glob

# Define the author names which did absolute measurments that are trustworth.
authors = ['schmidt2016', 'li2014', 'valgepea2013', 'peebo2015']

# Load inthe datasets and concatenate.
data = pd.concat([
    pd.read_csv(f'../../../data/{auth}_longform_annotated.csv') for auth in authors],
    sort=False)
data.drop(columns=['reported_tot_per_cell', 'reported_fg_per_cell', 
                   'annotation', 'Unnamed: 0', 'reported_volume', 
                   'corrected_volume'], axis=1, inplace=True)

data.to_csv('../../../data/compiled_absolute_measurements.csv', index=False)


# %%
