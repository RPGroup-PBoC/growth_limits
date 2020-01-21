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
bokeh.io.output_notebook()
color, palette = prot.viz.bokeh_theme()

# Load in the data
raw_data = pd.read_csv('../../data/schmidt2016_longform.csv')
sectors = pd.read_csv('../../data/schmidt2016_cog_class_sectoring.csv')


#%%
# Define the positions 
groups = sectors['group'].unique()

layout = []
for g1 in groups:
    row = []
    for g2 in groups:
        g1_sectors = sectors[sectors['group']==g1]
        g2_sectors = sectors[sectors['group']==g2]
        g1_sectors.sort_values('growth_rate_hr', inplace=True)
        g2_sectors.sort_values('growth_rate_hr', inplace=True)
        source = ColumnDataSource(
            {'g1_frac':g1_sectors['frac_mass'].values,
            'g2_frac':g2_sectors['frac_mass'].values,
            'condition':g1_sectors['condition']})
        p = bokeh.plotting.figure(width=200, height=200, tools=['hover'],
                    tooltips=[('condition', '@condition')])
        p.circle(x='g2_frac', y='g1_frac', source=source)
        row.append(p)
    layout.append(row)

# Format the axes labels:
for i, g in enumerate(groups):
    layout[i][0].yaxis.axis_label = g 
    layout[-1][i].xaxis.axis_label = g 

grid = bokeh.layouts.gridplot(layout)
bokeh.io.show(grid)

# %%
