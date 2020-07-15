#%%
import numpy as np
import pandas as pd 
import bokeh.io
from bokeh.models import *
import bokeh.layouts
import bokeh.plotting
import prot.viz
import tqdm
colors, palette = prot.viz.bokeh_theme()
dataset_colors = prot.viz.dataset_colors()

# ##############################################################################
# DATA CLEANING AND AGGREGATION
# ##############################################################################
# Load the three datasets
prots = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
cplx = pd.read_csv('../../../data/compiled_annotated_complexes.csv')
cplx = cplx[['gene_name', 'b_number', 'condition', 'growth_rate_hr', 'go_terms',
             'complex', 'complex_annotation', 'dataset', 'dataset_name',
             'n_subunits', 'n_units', 'gene_product']]

# Condense the complex data frame to an easily paresable form. 
cplx_numeric_dfs = []
cplx_desc_dfs = []

for g, d in tqdm.tqdm(cplx.groupby(['complex_annotation', 
                                    'complex', 'go_terms']), 
                                    desc='Condensing complex data sets...'):
    # Assemble the stoichiometric formula
    subunits, number = d['gene_name'].unique(), d['n_subunits'].unique()
    formula = ''
    for s, n in zip(subunits, number):
        formula += f'[{s[0].upper()}{s[1:]}]<sub>{int(n)}</sub>'

    # Assemble a descriptive data frame
    cplx_desc = pd.DataFrame([])
    cplx_desc = cplx_desc.append({'complex_annotation':g[0],
                                  'complex': g[1],
                                  'go_terms': g[2],
                                  'formula': formula}, 
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
                          'condition': _g[2],
                          'growth_rate_hr':_g[3],
                          'color': dataset_colors[_g[0]]
                          }, ignore_index=True)
    cplx_numeric_dfs.append(cplx_numeric)


# Define the column data sources
cplx_desc = pd.concat(cplx_desc_dfs, sort=False)
cplx_numeric = pd.concat(cplx_numeric_dfs, sort=False)

#%%
cplx_desc_source = ColumnDataSource(cplx_desc)
cplx_numeric_source = ColumnDataSource(cplx_numeric)


# ##############################################################################
# FIGURE CANVAS DECLARATION AND INTERACTION INSTANTIATION
# ##############################################################################
bokeh.io.output_file('./explorer.html')

# Define the menu for selecting complexes 
complex_menu = [(d, g) for g, d in zip(cplx['complex_annotation'].unique(), cplx['complex'].unique())]
complex_selection = Select(title='EcoCyc Complex Annotation', 
                           value='select annotation',
                           options=complex_menu)

# Define th menu for selecting proteins
prots.sort_values(by='gene_name', inplace=True)
protein_menu = [(d, g) for g, d in zip(prots['gene_name'].unique(), prots['b_number'].unique())]
protein_selection = Select(title='protein name', options=protein_menu,
                           value='select gene product')

# Define the selector for log or linear y scaling. 
scaling = RadioButtonGroup(labels=['log\u2081\u2080', 'linear'])

# Define the slector for min, max, or median complex abundance
agg_fn = RadioButtonGroup(labels=['minimum', 'maximum', 'median', 'mean'])

# Define the figure canvases
complex_canvas = bokeh.plotting.figure(width=400, height=400, 
                                       x_axis_label='growth rate [per hr]',
                                       y_axis_label='complex abundance per cell')
GO_canvas = bokeh.plotting.figure(width=400, height=400, 
                                  x_axis_label='growth rate [per hr]',
                                  y_axis_label='process member abundance per cell')
prot_canvas = bokeh.plotting.figure(width=400, height=400, 
                                  x_axis_label='growth rate [per hr]',
                                  y_axis_label='protein abundance per cell')


# Define description divs and tables. 
complex_description_field = Div(text="")
# ##############################################################################
# CALLBACK DEFINTITION AND ASSIGNMENT
# ##############################################################################

# Define the argument dictionaries
cplx_args = {'cplx_desc_source': cplx_desc_source, 
             'cplx_numeric_source':cplx_numeric_source, 
             'cplx_select_input': complex_selection,
             'scaling_input': scaling,
             'agg_method_input': agg_fn,
             'cplx_desc_field':complex_description_field
             }

cplx_cb = prot.viz.load_js('explorer_complex.js', args=cplx_args)
complex_selection.js_on_change('value', cplx_cb)

# Define the widget boxes 
complex_box = bokeh.layouts.row(scaling, agg_fn) 
complex_box = bokeh.layouts.column(complex_box, complex_selection)
go_box = bokeh.layouts.row(scaling, agg_fn)
prot_box = bokeh.layouts.column(scaling, protein_selection)
complex_layout = bokeh.layouts.row(complex_canvas, complex_box)
go_layout = bokeh.layouts.row(GO_canvas, go_box)
protein_layout = bokeh.layouts.row(prot_canvas, prot_box)

# Define the tabs
complex_tab = Panel(child=complex_layout, title='by functional complex')
GO_tab = Panel(child=go_layout, title='by GO annotation')
prot_tab = Panel(child=protein_layout, title='by protein')
tabs = Tabs(tabs=[complex_tab, GO_tab, prot_tab])
bokeh.io.save(tabs)

# %%
