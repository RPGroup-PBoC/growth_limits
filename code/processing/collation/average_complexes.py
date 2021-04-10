#%%
import numpy as np 
import pandas as pd
import tqdm 
data = pd.read_csv('../../../data/compiled_annotated_complexes.csv')

complex_df = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby(['condition', 'growth_rate_hr', 
                          'dataset', 'dataset_name', 'complex', 
                          'complex_annotation'])):
    if len(d) == 1:
        if d['n_subunits'].values[0] == 1:
            complex_type = 'monomeric' 
        else:
            complex_type = 'multimeric' 
    else:
       complex_type = 'complex' 

    proteins = list(d['gene_name'].values)
    subunits = list(d['n_subunits'].values.astype(int))
    formula = ''.join([f'[{prot}]_{sub}' for prot, sub in zip(proteins, subunits)]) 
    # Compute the mass in kilodaltons
    complex_mass = 602213665.167516E-3 * np.sum((d['fg_per_cell'].values / d['tot_per_cell'].values) * d['n_subunits'].values)
    min_number = d['n_units'].min()
    max_number = d['n_units'].max()
    mean_number = d['n_units'].mean()
    median_number = d['n_units'].median()
    cog_letter = d['cog_letter'].unique()[0]
    cog_class = d['cog_class'].unique()[0]
    cog_category= d['cog_category'].unique()[0]
    complex_df = complex_df.append({
                    'complex': g[-2],
                    'complex_annotation': g[-1],
                    'formula': formula,
                    'complex_mass_kda': complex_mass,
                    'cog_category': cog_category,
                    'cog_class': cog_class,
                    'cog_letter': cog_letter,
                    'dataset': g[2],
                    'dataset_name': g[3],
                    'condition': g[0],
                    'growth_rate_hr': g[1],
                    'strain': d['strain'].values[0],
                    'min_number': min_number,
                    'max_number': max_number,
                    'mean_number': mean_number,
                    'median_number': median_number,

    }, ignore_index=True) 


# %%
complex_df.to_csv('../../../data/compiled_complex_abundances.csv', index=False)
# %%
