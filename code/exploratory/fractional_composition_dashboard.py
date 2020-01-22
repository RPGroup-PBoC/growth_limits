#%%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models import * 
import prot.stats
import prot.viz
import bokeh.palettes
colors, palette = prot.viz.bokeh_theme()
np.random.seed(666)

# Define a function that will assign rectangular bounds / colors for treemaps
def group_and_assign(df, groupby='condition', key='frac_mass'):
    """
    Function to assign rectangular bounds for treemaps and colors each sector
    """
    dfs = []
    for g, d in df.groupby(groupby):
        _df = prot.viz.assign_rect_bounds(d, key)
        dfs.append(_df)
    return pd.concat(dfs)

def assign_color(df, key, palette):
    if len(df[key].unique()) > 200:
        replace = True
        n = 200
    else:
        replace = False
        n = len(df[key].unique()) + 2
    _palette = np.random.choice(palette(n), replace=replace, size=len(df[key].unique()))
    
    for i, k in enumerate(df[key].unique()):
        df.loc[df[key]==k, 'color'] = _palette[i]
    return df

def process(df, palette, groupby='condition', key='frac_mass', color_key='group'):
    bounded = group_and_assign(df, groupby=groupby, key=key)
    colored = assign_color(bounded, key=color_key, palette=palette)
    return colored

# %%
bokeh_palette = bokeh.palettes.viridis

# Load the COG sectoring data. 
class_sectors = pd.read_csv('../../data/schmidt2016_cog_class_sectoring.csv')
class_sectors = process(class_sectors, bokeh_palette)
class_sectors['display_color'] = class_sectors['color']
class_sectors['display_alpha'] = 0.25
desc_sectors = pd.read_csv('../../data/schmidt2016_cog_desc_sectoring.csv')
desc_sectors = process(desc_sectors, bokeh_palette, groupby=['condition', 'cog_class'])
desc_sectors['display_color'] = desc_sectors['color']
desc_sectors['display_alpha'] = 0.25
gene_sectors = pd.read_csv('../../data/schmidt2016_cog_gene_sectoring.csv')
gene_sectors = process(gene_sectors, bokeh_palette, groupby=['condition', 'cog_desc'])
gene_sectors['display_color'] = gene_sectors['color']
gene_sectors['display_alpha'] = 0.25 

#%%
bokeh.io.output_file('../../figs/fractional_composition_dashboard.html')

# Define the sectors as ColumnDataSources
class_tree_source = ColumnDataSource(class_sectors)
class_rate_source = ColumnDataSource(class_sectors)
desc_tree_source = ColumnDataSource(desc_sectors)
desc_rate_source = ColumnDataSource(desc_sectors)
gene_tree_source = ColumnDataSource(gene_sectors)
gene_rate_source = ColumnDataSource(gene_sectors)

# Define the condition selector
class_sectors.sort_values('growth_rate_hr', inplace=True)
condition_selector = Select(title="growth condition", options=list(class_sectors['condition'].unique()),
                            value='lb_miller')

# Define the condition, class, and subgroup filters 
condition_filter = GroupFilter(column_name='condition', group='lb_miller')
class_filter = GroupFilter(column_name='cog_class', group='')
desc_filter = GroupFilter(column_name='cog_desc', group='')


# Define the views
condition_view = CDSView(source=class_tree_source, filters=[condition_filter])
class_view = CDSView(source=desc_tree_source, filters=[condition_filter, 
                                                        class_filter])
gene_view = CDSView(source=gene_tree_source, filters=[condition_filter, 
                                                      class_filter, desc_filter])
class_rate_view = CDSView(source=desc_tree_source, filters=[class_filter])
gene_rate_view = CDSView(source=gene_tree_source, filters=[class_filter, desc_filter])


# Define the callbacks
class_click_cb = """
    var cog_class_ind = class_tree_source.selected['1d'].indices[0];
    var cog_class = class_tree_source.data['group'][cog_class_ind];
    class_filter.group = cog_class;
    class_view.filters = [condition_filter, class_filter];
    class_rate_view.filters = [class_filter];
    desc_tree_source.data.view = class_view
    desc_rate_source.data.view = class_rate_view;
    desc_tree_source.change.emit();

    // Update the color of the rate plot
    for (var i = 0; i < class_tree_source.data['group'].length; i++ ) {
        var sel_cog_class = class_tree_source.data['group'][i]
        if (sel_cog_class === cog_class) {
            class_rate_source.data['display_color'][i] = 'tomato';
            class_rate_source.data['display_alpha'][i] = 1;
            }
        else {
            class_rate_source.data['display_color'][i] = class_rate_source.data['color'][i];
            class_rate_source.data['display_alpha'][i] = 0.25;
        } 

    }
    class_rate_source.change.emit();
    desc_rate_source.change.emit();
    """

desc_click_cb = """
    var cog_desc_ind = desc_tree_source.selected['1d'].indices[0];
    var cog_desc = desc_tree_source.data['group'][cog_desc_ind];
    desc_filter.group = cog_desc;
    gene_view.filters = [condition_filter, class_filter, desc_filter];
    gene_rate_view.filters = [desc_filter];
    gene_tree_source.data.view = gene_view;
    gene_rate_source.data.view = gene_rate_view
    gene_tree_source.change.emit();
    gene_rate_source.change.emit();

    // Update the color of the rate plot
    for (var i = 0; i < desc_tree_source.data['group'].length; i++ ) {
        var sel_cog_desc = desc_tree_source.data['group'][i]
        if (sel_cog_desc === cog_desc) {
            desc_rate_source.data['display_color'][i] = 'tomato';
            desc_rate_source.data['display_alpha'][i] = 1;
            }
        else {
            desc_rate_source.data['display_color'][i] = desc_rate_source.data['color'][i];
            desc_rate_source.data['display_alpha'][i] = 0.25;
        } 

    }
    desc_rate_source.change.emit();
"""

gene_click_cb = """
    var gene_ind = gene_tree_source.selected['1d'].indices[0];
    var gene = gene_tree_source.data['group'][gene_ind]

    // Update the color of the rate plot
    for (var i = 0; i < gene_tree_source.data['group'].length; i++ ) {
        var sel_gene = gene_tree_source.data['group'][i]
        if (sel_gene === gene) {
            gene_rate_source.data['display_color'][i] = 'tomato';
            gene_rate_source.data['display_alpha'][i] = 1;
            }
        else {
            gene_rate_source.data['display_color'][i] = gene_rate_source.data['color'][i];
            gene_rate_source.data['display_alpha'][i] = 0.25;
        } 

    }
    gene_rate_source.change.emit();
 

"""

condition_select_cb = """
    condition_filter.group = selection.value;
    condition_view.filters = [condition_filter];
    class_tree_source.data.view = condition_view;
    class_tree_source.change.emit();

"""


# Assemble the callback
args = {'class_tree_source':class_tree_source,
        'class_rate_source':class_rate_source,
        'desc_tree_source':desc_tree_source,
        'desc_rate_source': desc_rate_source,
        'gene_tree_source':gene_tree_source,
        'gene_rate_source':gene_rate_source,
        'condition_filter':condition_filter,
        'class_filter':class_filter,
        'desc_filter':desc_filter,
        'condition_view':condition_view,
        'class_view':class_view,
        'class_rate_view':class_rate_view,
        'gene_rate_view':gene_rate_view,
        'gene_view':gene_view,
        'selection': condition_selector}


cb = CustomJS(args=args, 
                code=condition_select_cb + class_click_cb + desc_click_cb + gene_click_cb)


# Set up bokeh treemap axes
class_treemap = bokeh.plotting.figure(width=300, height=300, x_range=[0, 500],
                                      y_range=[0, 500], 
                                      title='proteome occupancy by COG class',
                                      tools=['wheel_zoom', 'pan', 'tap', 'hover'],
                                      tooltips=[('COG class', '@group')])

desc_treemap = bokeh.plotting.figure(width=300, height=300, x_range=[0, 500],
                                      y_range=[0, 500], 
                                      title='COG class occupancy by subgroup',
                                      tools=['wheel_zoom', 'pan', 'tap', 'hover'],
                                      tooltips=[('COG class', '@cog_class'),
                                                ('COG subgroup', '@group')])
gene_treemap = bokeh.plotting.figure(width=300, height=300, x_range=[0, 500],
                                      y_range=[0, 500], 
                                      title='COG subclass occupancy by gene',
                                      tools=['wheel_zoom', 'pan', 'tap',
                                      'hover'],
                                      tooltips=[('COG class', '@cog_class'),
                                                ('COG subgroup', '@cog_desc'),
                                                ('gene abbreviation', '@group'),
                                                ('gene name', '@desc')])
# Add the click callbacks
class_click_event = class_treemap.select(type=TapTool)
desc_click_event = desc_treemap.select(type=TapTool)
gene_click_event = gene_treemap.select(type=TapTool)
class_click_event.callback = cb
desc_click_event.callback = cb
gene_click_event.callback = cb

# Add the selection callback
condition_selector.js_on_change('value', cb)


# Set up the growth rate axes
class_rate = bokeh.plotting.figure(width=300, height=300,
                                title="total proteome occpancy (mass fraction)",
                                x_axis_label='growth rate [hr^-1]',
                                y_axis_label='mass fraction of proteome')
desc_rate = bokeh.plotting.figure(width=300, height=300,
                                title="COG class occupancy (mass fraction)",
                                x_axis_label='growth rate [hr^-1]',
                                y_axis_label='mass fraction of class')
gene_rate = bokeh.plotting.figure(width=300, height=300,
                                title="COG subgroup occupancy (mass fraction)",
                                x_axis_label='growth rate [hr^-1]',
                                y_axis_label='mass fraction of subgroup')

# Populate the treemaps with data. 
class_treemap.quad(bottom='bottom', top='top', left='left', right='right', 
                  source=class_tree_source, line_color='white', line_width=1, 
                  color='color', view=condition_view,
                  hover_color='tomato')
desc_treemap.quad(bottom='bottom', top='top', left='left', right='right', 
                  source=desc_tree_source, line_color='white', line_width=1, 
                  color='color', view=class_view,
                  hover_color='tomato')
gene_treemap.quad(bottom='bottom', top='top', left='left', right='right', 
                  source=gene_tree_source, line_color='white', line_width=1, 
                  color='color', view=gene_view,
                  hover_color='tomato')

# Populate the occupancy plots.
class_rate.circle(x='growth_rate_hr', y='frac_mass', color='display_color', 
                 alpha='display_alpha', source=class_rate_source, size=8, 
                 line_color='black') 
desc_rate.circle(x='growth_rate_hr', y='frac_mass', color='display_color', 
                alpha='display_alpha',
                source=desc_rate_source, view=class_rate_view, size=8,
                line_color='black', hover_color='tomato')
gene_rate.circle(x='growth_rate_hr', y='frac_mass', color='display_color', 
                alpha='display_alpha',
                source=gene_rate_source, view=gene_rate_view,
                size=8, line_color='black', hover_color='tomato')

# Add a DIV with instructions for use
text= """
This dashboard permits the growth-rate dependence of the various Clusters of Orthologous Groups (COG) classes and subgroups. It has multiple levels of interaction described below in order.

1) Growth Condition Selector. This dropdown allows you to select the particular growth condition. The entries are sorted in ascending order (slow-growing to fast).

2) Square "treemap" plots. The second row of plots shows the partitioning of the proteome into COG classes (left), COG subgroups (middle), and individual genes (right). The area of each square corresponds to the fraction of the section (by masss) which is occupied. Hovering over a section provides information about that subgroup. Clicking a square will break it further down into subgroup and genes.

3)
"""

# Define the layout object. 
treemap_row = bokeh.layouts.row(class_treemap, desc_treemap, gene_treemap)
rate_row = bokeh.layouts.row(class_rate, desc_rate, gene_rate)
plot = bokeh.layouts.column(condition_selector, treemap_row, rate_row)
bokeh.io.save(plot)
# %%


# %%
