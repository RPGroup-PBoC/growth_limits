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
color, palette = prot.viz.bokeh_theme()

# Define the color palettes
viridis = bokeh.palettes.viridis

# For each condition, compute the rect bounds and assign colors
def group_and_assign(df, groupby='condition', key='frac_mass', palette=viridis,
                    **kwargs):
    """
    Function to assign rectangular bounds for treemaps and colors each sector
    """
    dfs = []
    for g, d in df.groupby(groupby):
        _df = prot.viz.assign_rect_bounds(d, key, **kwargs)
        dfs.append(_df)
    return pd.concat(dfs)

def assign_color(df, key, palette):
    if len(df[key].unique()) > 254:
        _palette = np.random.choice(palette(256), replace=True, size=len(df[key].unique()))
    else:
        _palette = palette(len(df[key].unique()) + 2)
    for i, k in enumerate(df[key].unique()):
        df.loc[df[key]==k, 'color'] = _palette[i]
    return df

# Load the datasets
cog_class_df = pd.read_csv('../../data/schmidt2016_cog_class_sectoring.csv')
cog_class_df = group_and_assign(cog_class_df, groupby=['condition'])
cog_class_df = assign_color(cog_class_df, 'group', viridis)
cog_desc_df = pd.read_csv('../../data/schmidt2016_cog_desc_sectoring.csv')
cog_desc_df = group_and_assign(cog_desc_df, groupby=['condition', 'cog_class'])
cog_desc_df = assign_color(cog_desc_df, 'group', viridis)
cog_gene_df = pd.read_csv('../../data/schmidt2016_cog_gene_sectoring.csv')
cog_gene_df = assign_and_color(cog_gene_df, groupby=['condition', 'cog_desc'])
cog_gene_df = assign_color(cog_gene_df, 'group', viridis)


#%%
# Define the bokeh output
bokeh.io.output_file('./cog_dashboard.html')

# Define a selector for the conditions
cog_class_df.sort_values(by='growth_rate_hr', inplace=True)
selection = Select(title='growth condition', 
    value='lb_miller', options=list(cog_class_df['condition'].unique()))

# Define the data source. 
cog_class_source = ColumnDataSource(cog_class_df)
cog_desc_source = ColumnDataSource(cog_desc_df)
cog_gene_source = ColumnDataSource(cog_gene_df)

# Define the group filter on the condition. 
condition_filter = GroupFilter(column_name='condition', group='lb_miller')
cog_class_filter = GroupFilter(column_name='cog_class', group='')
cog_desc_filter = GroupFilter(column_name='cog_desc', group='')

# Define the view with the assigned filter
condition_view = CDSView(source=cog_class_source, filters=[condition_filter])
cog_desc_view = CDSView(source=cog_desc_source, filters=[condition_filter, 
                                                    cog_class_filter])
cog_gene_view = CDSView(source=cog_gene_source, filters=[condition_filter, 
                                                    cog_class_filter, 
                                                    cog_desc_filter])

tooltips = [('COG class', '@group'),
            ('Mass fraction', '@frac_mass'),
            ('Count fraction', '@frac_count'),
            ('Growth condition', '@condition'),
            ('Growth rate (hr^-1)', '@growth_rate_hr')]
class_treemap = bokeh.plotting.figure(width=400, height=400, 
               tools=['tap', 'hover', 'wheel_zoom', 'pan'], 
               tooltips=tooltips, title='Proteome occupancy by COG class',
               x_range=[0, 500], y_range=[0, 500])
desc_treemap = bokeh.plotting.figure(width=400, height=400,
               tools=['tap', 'hover', 'wheel_zoom', 'pan'], 
               tooltips = [('COG class', '@cog_class'),
                           ('COG subgroup', '@group'),
                           ('Mass fraction', '@frac_mass'),
                           ('Count fraction', '@frac_count'), 
                           ('Growth condition', '@condition'),
                           ('Growth rate (hr^-1)', '@growth_rate_hr')],
                title='COG class occupancy by subgroup',
                x_range=[0, 500], y_range=[0, 500])
gene_treemap = bokeh.plotting.figure(width=400, height=400,
               tools=['hover', 'wheel_zoom', 'pan'], 
               tooltips = [('COG class', '@cog_class'),
                           ('COG subgroup', '@cog_desc'),
                           ('Gene ID', '@group'),
                           ('Gene Name', '@desc'),
                           ('Mass fraction', '@frac_mass'),
                           ('Count fraction', '@frac_count'), 
                           ('Growth condition', '@condition'),
                           ('Growth rate (hr^-1)', '@growth_rate_hr')],
                title='Subgroup occupancy by protein gene product',
                x_range=[0, 500], y_range=[0, 500])

# Set the cog_class treemap
class_treemap.quad(top='top', bottom='bottom', left='left', right='right',
                source=cog_class_source, color='color', line_color='white', 
                view=condition_view)
desc_treemap.quad(top='top', bottom='bottom', left='left', right='right',
                  source=cog_desc_source, color='color', line_color='white',
                  view=cog_desc_view)
gene_treemap.quad(top='top', bottom='bottom', left='left', right='right',
                  source=cog_gene_source, color='color', line_color='white',
                  view=cog_gene_view)

for a in [class_treemap, desc_treemap, gene_treemap]:
    a.xaxis.major_tick_line_color = None
    a.yaxis.major_tick_line_color = None
    a.xaxis.major_label_text_font_size = "0pt"
    a.yaxis.major_label_text_font_size = "0pt"

# ##############################################################################

# ##############################################################################
# Interactivity
# ##############################################################################

# Identify the selection
class_click_cb = """
var cog_class_ind = cog_class_source.selected['1d'].indices[0];
var cog_class = cog_class_source.data['group'][cog_class_ind];
cog_class_filter.group = cog_class;
cog_desc_view.filters = [condition_filter, cog_class_filter];
cog_desc_source.data.view = cog_desc_view;
cog_desc_source.change.emit();
"""

desc_click_cb = """
var cog_desc_ind = cog_desc_source.selected['1d'].indices[0];
console.log(cog_desc_ind)
var cog_desc = cog_desc_source.data['group'][cog_desc_ind];
console.log(cog_desc)
cog_desc_filter.group = cog_desc;
cog_gene_view.filters = [condition_filter, cog_class_filter, cog_desc_filter];
cog_gene_source.data.view = cog_gene_view;
cog_gene_source.change.emit();
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
      'cog_gene_source':cog_gene_source, 'cog_desc_filter':cog_desc_filter,
      'cog_gene_view':cog_gene_view,
      'condition_filter':condition_filter, 'cog_class_filter':cog_class_filter, 
      'condition_view':condition_view, 'cog_desc_view':cog_desc_view,
      'selection':selection, 'cog_class_treemap':class_treemap}
selection_callback = CustomJS(args=args, code=selection_cb + class_click_cb + desc_click_cb)
class_click_callback = CustomJS(args=args, code=class_click_cb + desc_click_cb + selection_cb)
desc_click_callback = CustomJS(args=args, code=desc_click_cb + class_click_cb + selection_cb)

# Assign the callbacks to the interactions
class_click_event = class_treemap.select(type=TapTool)
desc_click_event = desc_treemap.select(type=TapTool)
class_click_event.callback = class_click_callback
desc_click_event.callback = desc_click_callback
selection.js_on_change('value', selection_callback)

row = bokeh.layouts.row(class_treemap, desc_treemap, gene_treemap)
col = bokeh.layouts.column(selection, row)
lay = bokeh.layouts.row(col)
bokeh.io.save(lay)

# %%
