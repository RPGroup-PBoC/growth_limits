#%%
import numpy as np
import pandas as pd
import networkx as nx
import bokeh.io
import bokeh.plotting
from bokeh.models import *
from bokeh.models.graphs import (from_networkx, NodesAndLinkedEdges, 
                                EdgesAndLinkedNodes, NodesOnly)
import prot.viz
import prot.stats
colors, palette = prot.viz.bokeh_theme()

# Load the complete gene-linked dataset
gene_uniprot_data = pd.read_csv('../../data/uniprot_biological_processes.csv')
gene_data = pd.read_csv('../../data/schmidt2016_longform.csv')

# Associate each gene and condition with a biological process
dfs = []
for g, d in gene_uniprot_data.groupby(['process', 'gene']):
   genes = gene_data[gene_data['gene']==g[1]].copy() 
   genes['biological_process'] = g[0]
   dfs.append(genes)
linked_genes = pd.concat(dfs)

# Define colors
process_colors = {p:np.random.choice(bokeh.palettes.d3['Category20b'][20]) for p in gene_uniprot_data['process'].unique()}

# Group by the condition and the process to compute the fractional mass
fractionated_processes = []
for g, d in linked_genes.groupby(['condition', 'growth_rate_hr']):
    # Get the total mass of the proteome in each condition, ignoring duplicates
    proteome_mass = gene_data[gene_data['condition']==g[0]]['fg_per_cell'].sum()

    # For each process, compute the fractional mass. **Note this wont sum to 1*
    process_mass = d.groupby('biological_process')['fg_per_cell'].sum().reset_index()
    process_mass['frac_mass'] = process_mass['fg_per_cell'].values / proteome_mass
    process_mass['condition'] = g[0]
    process_mass['growth_rate_hr'] = g[1]
    for p in d['biological_process'].unique():
        process_mass.loc[process_mass['biological_process']==p , 
                                          'highlight_color'] = process_colors[p]
    fractionated_processes.append(process_mass)
process_sectors = pd.concat(fractionated_processes)

#%%
# Pivot the dataframe to be untidy for simple correlation population
df = pd.DataFrame([])
for g, d in process_sectors.groupby('biological_process'):
    d.sort_values('growth_rate_hr', inplace=True)
    df[g] = d['frac_mass'].values
    
#%%
# Load the network dataset
network_data = pd.read_csv('../../data/uniprot_process_gene_network.csv')
network_data = network_data[network_data['n_genes'] != 0]

# Sort the network data by number of shared genes
network_data.sort_values('n_genes', inplace=True, ascending=False)
 
# Set up a list of nodes and edges. 
nodes = [(p1, {'process': p1, 
               'highlight_color':process_colors[p1],
               'display_color': 'grey'}) for p1 in gene_uniprot_data['process'].values]
edges = [(p1, p2, {'weight': j,
                   'width': 0.5,
                   'display_color': 'grey'}) for p1, p2, n in zip(network_data['process_1'].values,
                                       network_data['process_2'].values,
                                       network_data['n_genes'].values)]

# Convert the network data to a networkx object
graph = nx.Graph()
graph.add_nodes_from(nodes)
graph.add_edges_from(edges)

# %%
bokeh.io.output_file('../../figs/uniprot_process_gene_network.html')

# Define the graph
graph_renderer = from_networkx(graph, nx.spring_layout,  scale=25, center=(0, 0))

# Style the nodes
graph_renderer.node_renderer.glyph = Circle(size=8, fill_color='highlight_color', line_color='black', line_width=0.75)
graph_renderer.node_renderer.hover_glyph = Circle(fill_color=colors['red'],
                                                  line_color='black', 
                                                  line_width=0.75)
graph_renderer.node_renderer.selection_glyph = Circle(fill_color='highlight_color')

# Style the edges
graph_renderer.edge_renderer.glyph = MultiLine(line_color='display_color', 
                                               line_width=2)
graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=colors['red'], 
                                                     line_width='width')
graph_renderer.edge_renderer.selection_glyph = MultiLine(line_width='width')

# Set up the figure canvases. 
graph_plot = bokeh.plotting.figure(width=400, height=400, x_range=[-25, 25],
                        y_range=[-25, 25], 
                        tools=['hover', 'wheel_zoom', 'pan', 'tap'],
                        tooltips=[('UniProt Biological Process', '@process'),
                                  ('color', '$color[swatch]:highlight_color')])

# Set up the growth rate plot. 
correlation_plot = bokeh.plotting.figure(width=400, height=400, 
                                 x_axis_label='mass fraction of proteome (selected)',
                                 y_axis_label='mass fraction of proteome (compared)',
                                 tools=['hover', 'wheel_zoom', 'pan', 'tap'],
                                 tooltips=[('UniProt Biological Process', '@biological_process'),
                                           ('condition', '@condition'),
                                           ('growth rate [hr^-1]', '@growth_rate_hr'),
                                           ('mass fraction', '@frac_mass')])

# Set up the growth rate plot. 
rate_plot = bokeh.plotting.figure(width=800, height=400, 
                                 x_axis_label='growth rate [hr^-1]',
                                 y_axis_label='mass fraction of proteome',
                                 tools=['hover', 'wheel_zoom', 'pan', 'tap'],
                                 tooltips=[('UniProt Biological Process', '@biological_process'),
                                           ('condition', '@condition'),
                                           ('growth rate [hr^-1]', '@growth_rate_hr'),
                                           ('mass fraction', '@frac_mass')])


# Set up the data source and filters for the growth rate dependence
rate_source = ColumnDataSource(process_sectors)
rate_filter = IndexFilter(indices=[])
rate_view = CDSView(source=rate_source, filters=[rate_filter])

# Set up the data source for the correlation plot.
correlation_source =  {'selected':[], 'compared':[], 'color':[], 
                       'process':[], 'genes':[]}


# Populate the rate plot
rate_plot.circle(x='growth_rate_hr', y='frac_mass', 
                source=rate_source, color='highlight_color',
                size=8, line_color='black', line_width=0.75, view=rate_view)

# Populate the correlation plot
correlation_plot()

# Define a callback that changes colors upon selection
js = """
// Get the identity of the selected glyph
var ind = node_source.selected['1d'].indices[0];
var process_sel = node_source.data['process'][ind];
var process_color = node_source.data['highlight_color'][ind];

// Color the edge
var linked_nodes = [process_sel]
for (var i = 0; i < edge_source.data['start'].length; i++) {
    if ((edge_source.data['start'][i] === process_sel) || (edge_source.data['end'][i]===process_sel)) {
        var edge_start = edge_source.data['start'][i];
        var edge_end = edge_source.data['end'][i]; 
        
        // update the edge color
        edge_source.data['width'] = 2;

        // Change the edge color
        edge_source.data['display_color'][i] = process_color;

        // Keep track of what nodes are connected.
        if (edge_start == process_sel) { 
            linked_nodes.push(edge_end);
        }
        else {
            linked_nodes.push(edge_start);
        }

       }
    else {
        edge_source.data['display_color'][i] = 'grey';
    }
}

// Update the rate fraction plot. 
indices = []
for (var i = 0; i < rate_source.data['frac_mass'].length; i++) {
    if (linked_nodes.indexOf(rate_source.data['biological_process'][i]) >= 0) {
        indices.push(i);
    }
}

// Update the indices of the view
rate_filter.indices = indices;
rate_view.filters = [rate_filter];
rate_source.data.view = rate_view;


// Update the correlation plot
for (var i=0; i ; i++) {

}


// Emit changes to DOM
rate_source.change.emit();
edge_source.change.emit();
node_source.change.emit();
"""

args = {'node_source':graph_renderer.node_renderer.data_source,
        'edge_source':graph_renderer.edge_renderer.data_source,
        'rate_source':rate_source, 
        'rate_view':rate_view,
        'rate_filter':rate_filter}
cb = CustomJS(args=args, code=js)

# Add the callback to a click event
click_event = graph_plot.select(type=TapTool)
click_event.callback = cb

correlation = bokeh.plotting.figure() 
graph_plot.renderers.append(graph_renderer)
row = bokeh.layouts.row(graph_plot, correlation_plot)
layout = bokeh.layouts.column(row, rate_plot)
bokeh.io.save(layout)

# %%
