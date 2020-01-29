#%%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import scipy.spatial
bokeh.io.output_notebook()

# %%
# Load in the data of the proteome partitioning. 
cog_sectors = pd.read_csv('../../../data/schmidt2016_cog_class_sectoring.csv')

# %%
# Select a single condition
lb = cog_sectors[cog_sectors['condition']=='lb_miller']

# Randomly position the y coordinates for each class
y_pos = np.random.rand(len(lb['group'].unique()))
lb['y_pos'] = y_pos


# %%
vor = scipy.spatial.Voronoi(lb[['frac_mass', 'y_pos']].values)


# %%
# Set up a figure canvas. 
p = bokeh.plotting.figure(width=500, height=500)

p.circle(x='frac_mass', y='y_pos', source=lb)
# p.line(xvor.vertices)
bokeh.io.show(p)

# %%
vor.vertices

# %%
