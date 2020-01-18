#%%
import numpy as np
import pandas as pd
import networkx as nx
import bokeh.io
import bokeh.plotting
from bokeh.models import (ColumnDataSource, CDSView, GroupFilter, CustomJS, 
                        Circle, MultiLine, Plot)  
from bokeh.models.graphs import (from_networkx, NodesAndLinkedEdges, 
                                EdgesAndLinkedNodes, NodesOnly)
import prot.viz
colors, palette = prot.viz.bokeh_theme()

# Load the complete gene-linked dataset
gene_data = pd.read_csv('../../data/uniprot_biological_processes.csv')

# Load the network dataset
network_data = pd.read_csv('../../data/uniprot_process_gene_network.csv')
network_data = network_data[network_data['n_genes'] != 0]

# Set up a list of nodes and edges. 
nodes = [(p1, {'process': p1}) for p1 in gene_data['process'].values]
edges = [(p1, p2, {'weight':np.log10(n) + 0.1}) for p1, p2, n in zip(network_data['process_1'].values,
                                       network_data['process_2'].values,
                                       network_data['n_genes'].values)]

# Convert the network data to a networkx object

graph = nx.Graph()
graph.add_nodes_from(nodes)
graph.add_edges_from(edges)

# %%
bokeh.io.output_file('../../figs/uniprot_process_gene_network.html')

# Set up the figure canvas.
p = bokeh.plotting.figure(width=800, height=800, x_range=[-40, 40],
                        y_range=[-40, 40], 
                        tools=['hover', 'wheel_zoom', 'pan', 'tap'],
                        tooltips=[('UniProt Biological Process', '@process')])

# Define the graph
graph_renderer = from_networkx(graph, nx.spring_layout, scale=35, center=(0, 0), weight='weight')

# Style the nodes
graph_renderer.node_renderer.glyph = Circle(size=10, fill_color=colors['blue'], line_color='black', line_width=0.75)
graph_renderer.node_renderer.hover_glyph = Circle(size=14, fill_color=colors['green'], 
                                                  line_color='black', line_width=0.75)
graph_renderer.node_renderer.selection_glyph = Circle(size=14, fill_color=colors['red'], 
                                                  line_color='black', line_width=0.75)

# Style the edges
graph_renderer.edge_renderer.glyph = MultiLine(line_color='black')
graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=colors['green'], line_width=2)
graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=colors['red'], line_width=2)


graph_renderer.inspection_policy = NodesAndLinkedEdges()
graph_renderer.selection_policy = NodesAndLinkedEdges()
p.renderers.append(graph_renderer)
bokeh.io.save(p)

# %%
