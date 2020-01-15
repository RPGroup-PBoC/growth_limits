#%%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting 
import bokeh.layouts
import squarify as sq
import prot.viz
from bokeh.models import (ColumnDataSource, CDSView, GroupFilter, CustomJS,  
                          TapTool)
from bokeh.models.widgets import Dropdown, Select
import bokeh.palettes

# Define the color palettes
viridis = bokeh.palettes.viridis

# For each condition, compute the rect bounds and assign colors
def assign_and_color(df, groupby='condition', key='frac_mass', palette=viridis,
                    **kwargs):
    dfs = []
    for g, d in df.groupby(groupby):
        _df = prot.viz.assign_rect_bounds(d, key, **kwargs)
        if len(_df) > 254:
            _palette = np.random.choice(palette(256), replace=True, size=len(_df))
            _df['color'] = _palette
        else:
            _df['color'] = palette(len(_df) + 2)[:-2]

        dfs.append(_df)
    return pd.concat(dfs)

# Load the datasets
cog_class_df = pd.read_csv('../../data/schmidt2016_cog_class_sectoring.csv')
cog_class_df = assign_and_color(cog_class_df)
cog_desc_df = pd.read_csv('../../data/schmidt2016_cog_desc_sectoring.csv')
cog_desc_df = assign_and_color(cog_desc_df)
cog_gene_df = pd.read_csv('../../data/schmidt2016_cog_gene_sectoring.csv')
cog_gene_df = assign_and_color(cog_gene_df)

bokeh.io.output_file('./dashboard.html')
color, palette = prot.viz.bokeh_theme()

# Define a selector for the conditions
selection = Select(title='growth condition', 
    value='lb_miller', options=list(cog_class_df['condition'].unique()))

# Define the data source. 
cog_class_source = ColumnDataSource(cog_class_df)
cog_desc_source = ColumnDataSource(cog_desc_df)

# Define the group filter on the condition. 
condition_filter = GroupFilter(column_name='condition', group='lb_miller')
cog_class_filter = GroupFilter(column_name='cog_class', 
                group='INFORMATION STORAGE AND PROCESSING')

# Define the view with the assigned filter
condition_view = CDSView(source=cog_class_source, filters=[condition_filter])
cog_desc_view = CDSView(source=cog_desc_source, filters=[condition_filter, 
                                                    cog_class_filter])
tooltips = [('COG class', '@group'),
            ('Mass fraction', '@frac_mass'),
            ('Count fraction', '@frac_count')]
class_treemap = bokeh.plotting.figure(width=600, height=600, 
               tools=['tap', 'hover', 'wheel_zoom', 'pan'], 
               tooltips=tooltips)
desc_treemap = bokeh.plotting.figure(width=600, height=600,
               tools=['hover', 'wheel_zoom', 'pan'], 
               tooltips = [('COG class', '@cog_class'),
                           ('Description', '@group'),
                           ('Mass fraction', '@frac_mass'),
                           ('Count fraction', '@frac_count')])
# Set the cog_class treemap
class_treemap.quad(top='top', bottom='bottom', left='left', right='right',
                source=cog_class_source, color='color', line_color='white', 
                view=condition_view)
desc_treemap.quad(top='top', bottom='bottom', left='left', right='right',
                  source=cog_desc_source, color='color', line_color='white',
                  view=cog_desc_view)
# ##############################################################################
# Interactivity
# ##############################################################################

# Identify the selection
click_cb = """
var cog_class_ind = cog_class_source.selected['1d'].indices[0];
var cog_class = cog_desc_source.data['cog_class'][cog_class_ind];
console.log(cog_class)
cog_class_filter.group = cog_class;
cog_desc_view.filters = [condition_filter, class_filter];
cog_desc_source.data.view = cog_desc_view;
cog_desc_source.change.emit();
"""

# Add the selection
selection_cb = """
    // Update the main area plot upon changed grouping
    condition_filter.group = selection.value;
    condition_view.filters[0] = condition_filter;
    cog_class_source.data.view = condition_view;
    cog_class_source.change.emit();
"""
# Define a callback for updating the plot. 
args={'cog_class_source':cog_class_source, 'cog_desc_source':cog_desc_source, 
      'condition_filter':condition_filter, 'cog_class_filter':cog_class_filter, 
      'condition_view':condition_view, 'cog_desc_view':cog_desc_view,
      'selection':selection, 'cog_class_treemap':class_treemap}
selection_callback = CustomJS(args=args, code=selection_cb + click_cb)
click_callback = CustomJS(args=args, code=click_cb + selection_cb)

# Assign the callbacks to the interactions
click_event = class_treemap.select(type=TapTool)
click_event.callback = click_callback
selection.js_on_change('value', selection_callback)

col = bokeh.layouts.column(class_treemap, desc_treemap)
lay = bokeh.layouts.row(col, selection)
bokeh.io.save(lay)

# %%
