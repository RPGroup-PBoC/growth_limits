#%%
import numpy as np 
import pandas as pd 
import bokeh.io
import bokeh.plotting
from bokeh.models import * 
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.bokeh_theme()
bokeh.io.output_file('complex_protein_explorer.html')

# Load the three datasources
complex_desc_df = pd.read_csv('./complex_annotations.csv')
complexes_df  = pd.read_csv('./complexes_compressed.csv')

proteins_df = pd.read_csv('./proteins_compressed.csv')
proteins_df = proteins_df[proteins_df['protein']!='a']

# Define the data sources
desc = ColumnDataSource(complex_desc_df)
complexes = ColumnDataSource(complexes_df)
proteins = ColumnDataSource(proteins_df)

# Define the two display sources
complex_display = ColumnDataSource({'x':[], 'y':[], 'l':[], 'c':[], 'cond':[]})
protein_display = ColumnDataSource({'x':[], 'y':[], 'l':[], 'c':[], 'cond':[]})

# Define the interactions
COG_selector = Select(options=list(complex_desc_df['class'].unique()), 
                        value='Select COG functional class')
complex_selector = Select(options=[], value='Select complex')
protein_selector = Select(options=list(np.sort(proteins_df['protein'].unique())),
                          value='Select protein product')
agg_fn = RadioButtonGroup(labels=['minimum', 'maximum', 'median', 'mean'],
                                active=3)
TOOLTIPS = [('growth rate [per hr]', '@x'),
            ('abundance [per cell]', '@y{int}'),
            ('growth condition', '@cond'),
            ('reference', '@l')]
# Define simple annotations
COG_title = Div(text="<b>Clusters of Orthologous Groups (COG) Class</b>")
complex_title = Div(text="<b>EcoCyc Complex Annotation</b>")
agg_title = Div(text="<b>Complex Abundance Aggregation Function</b>")
protein_title = Div(text="<b>Gene Name</b>")
complex_desc = Div(text="")
protein_desc = Div(text="")

# Define the canvases
p_complex = bokeh.plotting.figure(width=500, height=400, 
                                  x_axis_label='growth rate [per hr]',
                                  y_axis_label='complex abundance [per cell]',
                                  y_axis_type='log', tooltips=TOOLTIPS)
p_protein = bokeh.plotting.Figure(width=500, height=400,  
                                  x_axis_label='growth rate [per hr]',
                                  y_axis_label='protein abundance [per cell]',
                                  y_axis_type='log', tooltips=TOOLTIPS)
for p in [p_complex, p_protein]:
    p.y_range.range_padding = 1 
    p.y_range.range_padding_units = 'percent'

# Populate the canvases
p_complex.circle(x='x', y='y', color='c', legend_field='l', 
                line_color='black', size=10, source=complex_display)
p_protein.circle(x='x', y='y', color='c', legend_field='l', 
                line_color='black', size=10, source=protein_display)

# Define the callbacks
cog_args = {'COG_selector':COG_selector,
           'complex_selector': complex_selector,
           'desc':desc}
cog_cb =  prot.viz.load_js('cog_selector.js', args=cog_args)
COG_selector.js_on_change('value', cog_cb)

complex_args = {'complex_display':complex_display,
               'complexes':complexes,
               'complex_selector':complex_selector,
               'agg_fn':agg_fn,
               'desc':desc,
               'complex_desc':complex_desc}
complex_cb = prot.viz.load_js('complex_selector.js', args=complex_args)
complex_selector.js_on_change('value', complex_cb)
agg_fn.js_on_change('active', complex_cb)

protein_args = {'protein_display':protein_display,
                'protein_selector': protein_selector,
                'proteins':proteins,
                'protein_desc':protein_desc}
protein_cb = prot.viz.load_js('./protein_selector.js', args=protein_args)
protein_selector.js_on_change('value', protein_cb)
complex_box = bokeh.layouts.column(COG_title, COG_selector, complex_title, complex_selector, 
                            agg_title, agg_fn, complex_desc)
prot_box = bokeh.layouts.column(protein_title, protein_selector, protein_desc)
complex_row = bokeh.layouts.row(p_complex, complex_box)
prot_row = bokeh.layouts.row(p_protein, prot_box)
complex_tab = Panel(child=complex_row, title='Explore by Complex Annotation')
prot_tab = Panel(child=prot_row, title='Explore by Individual Protein')
tabs = Tabs(tabs=[complex_tab, prot_tab])
bokeh.io.save(tabs)
# %%
