#%%
import numpy as np
import pandas as pd
import tqdm
import prot.viz
dataset_colors = prot.viz.dataset_colors()
# ##############################################################################
# DATA CLEANING AND AGGREGATION
# ##############################################################################

# Set up a dictionary to translate the conditions to something meaningful
condition_dict ={
        'lb_miller': 'LB Miller',
        'rich': 'EZ MOPS Rich Defined Medium',
        'glucose_minimal': 'M9 Minimal Medium + Glucose',
        'MOPS complete without methionine': 'MOPS Complete (- Met.)',
        'MOPS complete': 'MOPS Complete',
        'chemostat_u0.5': 'M9 Minimal Medium, Chemostat (µ=0.5)',
        'chemostat_u0.35': 'M9 Minimal Medium, Chemostat (µ=0.35)',
        'chemostat_u0.2': 'M9 Minimal Medium, Chemostat (µ=0.2)',
        'chemostat_u0.12': 'M9 Minimal Medium, Chemostat (µ=0.12)',
        'acetate': 'M9 Minimal Medium + Acetate',
        'fumarate': 'M9 Minimal Medium + Fumarate',
        'galactose': 'M9 Minimal Medium + Galactose',
        'glucose': 'M9 Minimal Medium + Glucose',
        'glucosamine': 'M9 Minimal Medium + Glucosamine',
        'glycerol': 'M9 Minimal Medium + Glycerol',
        'pyruvate': 'M9 Minimal Medium + Pyruvate',
        'succinate': 'M9 Minimal Medium + Succinate',
        'xylose': 'M9 Minimal Medium + Xylose',
        'glycerol_pAA': 'M9 Minimal Medium + Glycerol + Amino Acids',
        'mannose': 'M9 Minimal Medium + Mannose',
        'fructose': 'M9 Minimal Medium + Fructose',
        'pH6': 'M9 Minimal Medium + Glucose (pH = 6)',
        'MOPS minimal': 'MOPS Minimal Medium',
        '42C': 'M9 Minimal Medium + Glucose (42° C)',
        'osmotic_stress_glucose': 'M9 Minimal Medium + Glucose (Osmotic Stress)'}
       
        
# Load the three datasets
prots = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
cplx = pd.read_csv('../../../data/compiled_annotated_complexes.csv')
cplx = cplx[['gene_name', 'b_number', 'condition', 'growth_rate_hr', 'go_terms',
             'complex', 'complex_annotation', 'dataset', 'dataset_name',
             'n_subunits', 'n_units', 'gene_product']]
cplx = cplx[cplx['complex_annotation'] != 'none assigned']
cplx = cplx[cplx['complex'] != 'none assigned']
cplx

# Condense the complex data frame to an easily paresable form. 
cplx_numeric_dfs = []
cplx_desc_dfs = []

for g, d in tqdm.tqdm(cplx.groupby(['complex_annotation', 'complex']), 
                                    desc='Condensing complex data sets...'):
    cplx_desc = pd.DataFrame([])
    for _g, _d in d.groupby(['gene_name', 'n_subunits', 'gene_product']):
        _list = [g[0], _g[-1]]
        for i, k in enumerate(_list):
            k = k.replace('<sup>', '')
            k = k.replace('</sup>', '')
            k = k.replace('<sub>', '')
            k = k.replace('</sub>', '')
            k = k.replace('<SUP>', '')
            k = k.replace('</SUP>', '')
            k = k.replace('<SUB>', '')
            k = k.replace('</SUB>', '')

            # Remove itaicization and other spans
            k = k.replace('<i>', '')
            k = k.replace('</i>', '')
            k = k.replace('<I>', '')
            k = k.replace('</I>', '')
            k = k.replace('<small>', '')
            k = k.replace('</small>', '')

            # Replace greeks
            k = k.replace('&alpha;', 'α')
            k = k.replace('&beta;', 'β')
            k = k.replace('&gamma;', 'γ')
            k = k.replace('&epsilon;', 'ε')
            k = k.replace('&theta;', 'θ')
            k = k.replace('&Psi;', 'Ψ')
            k = k.replace('&psi;', 'ψ')
            k = k.replace('&chi;', 'χ')
            k = k.replace('&Chi;', 'Χ')

            # Replace others
            k = k.replace('&mdash;', '-')
            _list[i] = k

        # Assemble a descriptive data frame
        cplx_desc = cplx_desc.append({'complex_annotation':_list[0],
                                  'complex': g[1],
                                  'protein': _g[0][0].upper() + _g[0][1:],
                                  'subunits': _g[1], 
                                  'func': _list[1]},
                                  ignore_index=True)

    cplx_desc_dfs.append(cplx_desc)
    # Iterate through each complex, dataset, and condition, and compute the aggs
    cplx_numeric = pd.DataFrame([])
    for _g, _d in d.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
        cplx_numeric = cplx_numeric.append({
                          'min':_d['n_units'].min(),
                          'max':_d['n_units'].max(),
                          'mean':_d['n_units'].mean(),
                          'median':_d['n_units'].median(),
                          'complex': g[1],
                          'dataset': _g[0],
                          'dataset_name':_g[1],
                          'condition': condition_dict[_g[2]],
                          'cond':_g[2],
                          'growth_rate_hr':_g[3],
                          'color': dataset_colors[_g[0]]
                          }, ignore_index=True)
    cplx_numeric_dfs.append(cplx_numeric)


# Define the column data sources
cplx_desc = pd.concat(cplx_desc_dfs, sort=False)
cplx_numeric = pd.concat(cplx_numeric_dfs, sort=False)
cplx_desc.to_csv('./cplx_desc.csv', index=False)
cplx_numeric.to_csv('./cplx_numeric.csv', index=False)

# %%
