#%%
import pandas as pd
file = pd.read_csv('schmidt2016_dataset.csv')
rates = pd.read_csv('schmidt2016_rates.csv')
# %% 
# Generate sets of keys for renaming
fg_keys = [key for key in file.keys() if '_fg' in key]
cv_keys = [key for key in file.keys() if '_cv' in key]
tot_keys = [key for key in file.keys() if '_tot' in key]
others = [key for key in file.keys() if ('_fg' not in key) &
                                        ('_cv' not in key) &
                                        ('_tot' not in key)]
keys = {'fg_per_cell': fg_keys, 'coeff_var': cv_keys, 'tot_per_cell':tot_keys} 

# Define the renamed conditions.
renamed_cols = {'Uniprot Accession': 'uniprot', 'Description':'desc', 
               'Bnumber':'b_number', 'Gene':'gene', 
               'Annotated functional COG groups (letter)': 'cog_class_letter',
               'Annotated functional COG group (description)': 'cog_desc',
               'Annotated functional COG class': 'cog_class', 
               'Peptides.used.for.quantitation': 'peptide',
               'Confidence.score':'confidence_score', 
               'Molecular weight (Da)': 'mw_da', 
               'Dataset': 'dataset'}

renamed_conds = {'Glucose': 'glucose',
               'LB': 'lb_miller', 
               'Glycerol + AA': 'glycerol_pAA',
               'Pyruvate': 'pyruvate',
               'Chemostat µ=0.5': 'chemostat_u0.5',
               'Chemostat µ=0.35': 'chemostat_u0.35',
               'Chemostat µ=0.20': 'chemostat_u0.2',
               'Chemostat µ=0.12': 'chemostat_u0.12',
               'Osmotic-stress glucose': 'osmotic_stress_glucose',
               '42°C glucose':'glucose_42C',
               'pH6 glucose': 'glucose_pH6',
               'Acetate': 'acetate',
               'Xylose': 'xylose',
               'Mannose': 'mannose',
               'Galactose': 'galactose',
               'Succinate': 'succinate',
               'Fructose': 'fructose',
               'Fumarate': 'fumarate',
               'Glucosamine': 'glucosamine',
               'Glycerol': 'glycerol',
               'stationary_1day': 'stationary_1day',
               'stationary_3day': 'stationary_3day'}
# %% ITerate through each data set and make a new data frame with proper names. 
dfs = []
for quants, vals in keys.items():
    _df = pd.DataFrame([])
    for k in vals:
        _split = k.split('_')[-1]
        split = k.split(f'_{_split}')[0]
        _df[quants] = file[k]
        _df['condition'] = renamed_conds[split]
        for _k, _v, in renamed_cols.items():
            _df[_v] = file[_k]
        dfs.append(_df)
df = pd.concat(dfs)
df.to_csv('schmidt2016_longform.csv', index=False)
# %%
