#%%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting 
import bokeh.layouts
import squarify as sq
import prot.viz
import prot.stats
from bokeh.models import (ColumnDataSource, CDSView, GroupFilter, CustomJS,  
                          TapTool)
from bokeh.models.widgets import Dropdown, Select
import bokeh.palettes
colors, palette = prot.viz.bokeh_theme()

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
            _palette = np.random.choice(palette(len(_df) + 10), 
                                replace=True, size=len(_df))
            _df['color'] = _palette

        dfs.append(_df)
    return pd.concat(dfs)

# %%
# Load the coarse and fine-grained data set
processes = pd.read_csv('../../data/schmidt2016_uniprot_process_sectoring.csv')
genes = pd.read_csv('../../data/schmidt2016_uniprot_process_gene_sectoring.csv')

# Build the treemap for the genes
genes = assign_and_color(genes, groupby=['condition', 'uniprot_bio_process'])
#%%
# Instantiate the output file. 
bokeh.io.output_file('../../figs/biological_process_dashboard.html')

# Define the data sources
growth_source = ColumnDataSource(processes)
treemap_source = ColumnDataSource(genes)

# Define a view filter based on a selection
selection = Select(title="Biological Process (UniProt):", value='',
                   options=list(processes['uniprot_bio_process'].unique()))


# Define a view and group filter for the selected process
process_filter = GroupFilter(column_name='uniprot_bio_process', group='')
condition_filter = GroupFilter(column_name='condition', group='')
growth_view = CDSView(source=growth_source, filters=[process_filter])
treemap_view = CDSView(source=treemap_source, filters=[process_filter, condition_filter])

growth_tooltips = [('growth condition', '@condition'),
                    ('growth rate [hr^-1]', '@growth_rate_hr')]
treemap_tooltips = [('gene', '@group'),
                    ('growth condition', '@condition'),
                    ('COG class', '@cog_class'),
                    ('Description', '@desc')]
# Set up the canvases
treemap = bokeh.plotting.figure(width=450, height=450, x_range=[0, 500],
                                y_range=[0, 500], 
                                title='OCCUPANCY OF PROCESS BY PROTEIN',
                                tools=['wheel_zoom', 'pan','hover'],
                                tooltips=treemap_tooltips)
growth_plot = bokeh.plotting.figure(width=450, height=450, 
        x_axis_label='growth rate [hr^-1]', y_axis_label='mass fraction of proteome',
        x_range = [-0.05, 2], y_axis_type='log',
        tools=['wheel_zoom','pan', 'tap', 'hover'], tooltips=growth_tooltips)

# Populate the axes.
treemap.quad(bottom='bottom', left='left', top='top', right='right', color='color',
            source=treemap_source, view=treemap_view, line_color='white',
            alpha=0.85)
growth_plot.circle(x='growth_rate_hr', y='frac_mass', source=growth_source, color=colors['red'],
                    view=growth_view, size=7, alpha=0.75)

# Hide the ticks on the treemap
treemap.xaxis.major_tick_line_color = None
treemap.yaxis.major_tick_line_color = None
treemap.xaxis.major_label_text_font_size = "0pt"
treemap.yaxis.major_label_text_font_size = "0pt"
# Define the callback arguments
args = {
    'selection':selection,
    'process_filter':process_filter,
    'condition_filter':condition_filter,
    'growth_view': growth_view, 'treemap_view':treemap_view,
    'growth_source':growth_source, 'treemap_source':treemap_source}

# Define the callback codes
selection_code = """
    process_filter.group = selection.value;
    growth_view.filters[0] = process_filter;
    growth_source.data.view = growth_view;
    growth_source.change.emit();
    """
click_code = """
    var condition_ind = growth_source.selected['1d'].indices[0];
    console.log(condition_ind);
    var condition = growth_source.data['condition'][condition_ind]
    console.log(condition);
    process_filter.group = selection.value;
    condition_filter.group = condition;
    treemap_view.filters = [process_filter, condition_filter];
    treemap_source.data.view = treemap_view;
    treemap_source.change.emit();
"""

# Define the callbacks
selection_cb = CustomJS(args=args, code=selection_code)
click_cb = CustomJS(args=args, code=click_code)

# Assign callbacks as necessary
selection.js_on_change('value', selection_cb)
click_event = growth_plot.select(type=TapTool)
click_event.callback = click_cb

# Set up the layout
row = bokeh.layouts.row(growth_plot, treemap)
col = bokeh.layouts.column(selection, row)
lay = col
bokeh.io.save(lay)
 
# %%
