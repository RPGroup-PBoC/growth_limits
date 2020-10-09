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
# Load the data
cplx_prot_numeric_df = pd.read_csv('./cplx_prot_numeric.csv')
cplx_prot_desc_df = pd.read_csv('./cplx_prot_desc.csv')
cplx_prot_numeric_source = ColumnDataSource(cplx_prot_numeric_df)
cplx_prot_desc_source = ColumnDataSource(cplx_prot_desc_df)
prot_numeric_df = pd.read_csv('./prot_numeric.csv')
prot_desc_df = pd.read_csv('prot_desc.csv')
prot_numeric_source = ColumnDataSource(prot_numeric_df)
prot_desc_source = ColumnDataSource(prot_desc_df)

# Define the data sources
cplx_df = pd.read_csv('../../../data/compiled_annotated_complexes.csv')
cplx_desc_df = pd.read_csv('./cplx_desc.csv')
cplx_numeric_df = pd.read_csv('./cplx_numeric.csv')
cplx_source_data = ColumnDataSource(cplx_df)
cplx_desc_source = ColumnDataSource(cplx_desc_df)
cplx_numeric_source = ColumnDataSource(cplx_numeric_df)

# Define the display sources. 
cplx_display_source = ColumnDataSource({'x':[], 'y':[], 'c':[], 'l':[],
                                        'condition': [], 'cond':[]})
cplx_table_source = ColumnDataSource({'protein':[], 'subunits':[], 
                                      'relative_subunits':[],
                                      'observed':[],
                                      'func':[]})
cplx_prot_display_source = ColumnDataSource({'x':[], 'y':[], 'c':[], 'l':[],
                                         'condition':[]})
prot_display_source = ColumnDataSource({'x':[], 'y':[], 'c':[], 'l':[],
                                         'condition':[]})


# ##############################################################################
# FIGURE CANVAS DECLARATION AND INTERACTION INSTANTIATION
# ##############################################################################
bokeh.io.output_file('./data_explorer.html')

# Define the menu for selecting complexes 
cog_menu = [g for g in np.sort(cplx_desc_df['cog'].unique())]
complex_menu = {d:g for g, d in zip(cplx_desc_df['complex_annotation'].values, cplx_desc_df['complex'].values)}
complex_menu = [(k, v) for k, v in complex_menu.items()]
class_selection = Select(value='select COG functional class',
                         options=cog_menu,
                         width=400)
complex_selection = Select(
                           value='select annotation',
                           options=[],
                           width=400)

# Define the menu for selecting proteins
protein_selection = Select(options=[], value='select gene product')
protein_name_selection = Select(options=[], value='select gene product')

# Define the slector for min, max, or median complex abundance
agg_fn = RadioButtonGroup(labels=['minimum', 'maximum', 'median', 'mean'],
                         active=3)


# Define the hover for callbacks and interactions.
TOOLTIPS = [('growth rate [per hr]', '@x'),
            ('abundance [per cell]', '@y{int}'),
            ('growth condition', '@condition'),
            ('data source', '@l')]
            
# Define the figure canvases
complex_canvas = bokeh.plotting.figure(width=500, height=400, 
                                       x_axis_label='growth rate [per hr]',
                                       y_axis_label='complex abundance per cell',
                                       y_axis_type='log',
                                       tooltips=TOOLTIPS)
complex_canvas.y_range.range_padding = 1 
complex_canvas.y_range.range_padding_units = 'percent'
prot_canvas = bokeh.plotting.figure(width=500, height=400, 
                                  x_axis_label='growth rate [per hr]',
                                  y_axis_label='protein abundance per cell',
                                  y_axis_type='log', 
                                  tooltips=TOOLTIPS)
prot_canvas.y_range.range_padding = 1 
prot_canvas.y_range.range_padding_units = 'percent'

# Populate the canvases.
complex_canvas.circle(x='x', y='y', color='c', legend_field='l',
                     source=cplx_display_source, line_color='black', size=10)
prot_canvas.circle(x='x', y='y', color='c', legend_field='l',
                   source=cplx_prot_display_source, line_color='black', size=10)

# Define description divs and tables. 
class_title = Div(text='<b>Clusters of Orthologous Groups (COG) class</b>')
selector_title = Div(text='<b> EcoCyc complex annotation</b>')
agg_title = Div(text='<b>complex abundance aggregation method</b>')
prot_title = Div(text='<b>EcoCyc primary gene name</b>')
complex_description_field = Div(text="")
protein_description_field = Div(text="")
complex_table_cols =  [
    TableColumn(field='protein', title='protein name'),
    TableColumn(field='relative_subunits', title='expected stoichiometry'),
    TableColumn(field='observed', title='observed stoichiometry')
]
complex_table = DataTable(columns=complex_table_cols, source=cplx_table_source,
                          width=450)
# ##############################################################################
# CALLBACK DEFINTITION AND ASSIGNMENT
# ##############################################################################

# Arguments for COG class selection widget
class_args = {'class_select_input': class_selection,
              'complex_select_input': complex_selection,
              'cplx_desc_source': cplx_desc_source,
              'prot_desc_source': cplx_prot_desc_source,
              'prot_select_input':protein_selection}

# Arguments for complex selection
cplx_args = {'cplx_desc_source': cplx_desc_source, 
             'cplx_numeric_source':cplx_numeric_source, 
             'cplx_select_input': complex_selection,
             'cplx_table_source': cplx_table_source,
             'agg_method_input': agg_fn,
             'cplx_desc_field':complex_description_field,
             'cplx_display_source': cplx_display_source}

# Arguments for hover callback
cplx_hover_args = { 'cplx_desc_source': cplx_desc_source,
                   'cplx_source_data':cplx_source_data, 
                   'cplx_select_input':complex_selection,
                   'cplx_table_source': cplx_table_source,
                   'cplx_display_source':cplx_display_source}

# Arguments for protein selection.
cplx_prot_args = {'prot_numeric_source': cplx_prot_numeric_source,
             'prot_display_source': cplx_prot_display_source,
             'prot_select_input': protein_selection,
             'prot_div': protein_description_field}
prot_args = {'prot_numeric_source': prot_numeric_source,
             'prot_display_source': prot_display_source,
             'prot_select_input': protein_selection,
             'prot_div': protein_description_field}

class_cb = prot.viz.load_js('class_selection.js', args=class_args)
hover_cb = prot.viz.load_js('explorer_subunits.js', args=cplx_hover_args)
cplx_cb = prot.viz.load_js('explorer_complex.js',args=cplx_args)
prot_cb = prot.viz.load_js('prot_explorer.js', args=cplx_prot_args)
prot_name_cb = prot.viz.load_js('prot_name_explorer.js', args=prot_args)
class_selection.js_on_change('value', class_cb)
protein_selection.js_on_change('value', prot_cb)
complex_selection.js_on_change('value', cplx_cb)
agg_fn.js_on_change('active', cplx_cb)
complex_canvas.hover.callback = hover_cb

# Define the widget boxes 
complex_box = bokeh.layouts.column(class_title, class_selection, selector_title,
                                   complex_selection, agg_title, agg_fn,
                                   complex_description_field, complex_table)

prot_box = bokeh.layouts.column(class_title, class_selection, prot_title, protein_selection,
                                protein_description_field)
complex_plots = bokeh.layouts.column(complex_canvas)
complex_layout = bokeh.layouts.row(complex_plots, complex_box)
protein_layout = bokeh.layouts.row(prot_canvas, prot_box)

# Define the tabs
complex_tab = Panel(child=complex_layout, title='Explore By Complex Annotation')
prot_tab = Panel(child=protein_layout, title='Explore By Individual Protein')
tabs = Tabs(tabs=[complex_tab, prot_tab])
bokeh.io.save(tabs)

# %%
