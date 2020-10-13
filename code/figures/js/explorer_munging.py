#%%
import numpy as np 
import pandas as pd
import numpy as np 
import prot.viz
import tqdm
dataset_colors = prot.viz.dataset_colors()
# Define necessary renaming function and dict 
def rename(k):
    k = k.replace('<sup>', '')
    k = k.replace('</sup>', '')
    k = k.replace('<sub>', '')
    k = k.replace('</sub>', '')
    k = k.replace('<SUP>', '')
    k = k.replace('</SUP>', '')
    k = k.replace('<SUB>', '')
    k = k.replace('</SUB>', '')

    # Remove italicization and other spans
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
    return k
# Set up a dictionary to translate the conditions to something meaningful
condition_dict ={
        'lb_miller': 'LB Miller',
        'rich': 'EZ MOPS Rich Defined Medium',
        'glucose_minimal': 'Glucose Minimal Medium.',
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
#  COG Letter Dict
cog_dict = {
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Post-translational modification, protein turnover, and chaperones',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense Mechanisms',
    'W': 'Extracellular structures',
    'Z': 'Cytoskeleton',
    'A': 'RNA processsing and modification',
    'B': 'Chromatin structure and dynamics',
    'J': 'Translation, ribosomal structure, and biogenesis',
    'K': 'Transcription', 
    'L': 'Replication, recombination, and repair',
    'C': 'Energy production and conversion',
    'E': 'Aminio acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism', 
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism', 
    'P': 'Inorganic ion transport and metabolism', 
    'Q': 'Secondary metabolites, biosynthesis, transport, and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'X': 'Not Assigned',
    'Not Assigned': 'No Specific COG Annotation'}
        
# Load the individual data sets
complexes = pd.read_csv('../../../data/compiled_annotated_complexes.csv')
complexes = complexes[complexes['complex'] !='none assigned']
proteins = pd.read_csv('../../../data/compiled_absolute_measurements.csv')

# Set up a simplified list of complexes and descriptions
cplx_desc = pd.DataFrame([])
for g, d in tqdm.tqdm(complexes.groupby(['cog_letter', 'complex', 'complex_annotation'])):
    cplx_desc = cplx_desc.append({'class': cog_dict[g[0]], 'complex':g[1],
                                  'complex_annotation': rename(g[-1])},
                                  ignore_index=True)
cplx_desc.to_csv('./complex_annotations.csv', index=False)

#%%  Condense the complex data
grouped = complexes.groupby(['complex_annotation', 'dataset', 'dataset_name', 'condition', 'growth_rate_hr']
                             )['n_units'].agg(('min', 'max', 'median', 'mean')).reset_index()
grouped['color'] = [dataset_colors[k] for k in grouped['dataset'].values]
grouped['condition'] = [condition_dict[k] for k in grouped['condition'].values]
grouped['complex_annotation'] = [rename(str(k)) for k in grouped['complex_annotation'].values]
grouped.to_csv('./complexes_compressed.csv', index=False)

#%%
# Make a compressed dataframe of the protein abundances
products = [rename(str(k)) for k in proteins['gene_product'].values]
_colors = [dataset_colors[k] for k in proteins['dataset'].values]
_conditions = [condition_dict[k] for k in proteins['condition'].values]
prots_numeric = pd.DataFrame([])
prots_numeric['protein'] = proteins['gene_name']
prots_numeric['growth_rate_hr'] = proteins['growth_rate_hr']
prots_numeric['gene_product'] = products 
prots_numeric['condition'] = _conditions
prots_numeric['abundance'] = proteins['tot_per_cell']
prots_numeric['dataset_name'] = proteins['dataset_name']
prots_numeric['color'] = _colors
prots_numeric.to_csv('./proteins_compressed.csv', index=False)



# %%