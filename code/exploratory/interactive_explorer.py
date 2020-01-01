#%%
#%% [markdown
# # Exploring the Heinneman Data
#
# I've gone through and conerted some of the Heinneman data into a tidy-format
# to make exploration a bit easier. The point of this notebook is to build a toy
# interactive explorer using the Panel dashboarding framework to make things a
# bit simpler. 
# 
# %%
import numpy as np 
import pandas as pd 
import bokeh.io 
import bokeh.plotting 
import holoviews as hv
import panel as pn
pn.extension()
bokeh.io.output_notebook()
hv.extension('bokeh')

# %%
# Load the data 
data = pd.read_csv('../../data/schmidt2016_longform.csv')

#%%
# %opts Scatter [width = 600, height = 400]
hv.Scatter(data, kdims=['growth_rate_hr'], vdims=['tot_per_cell',
                'cog_desc']).groupby(['cog_desc'])



# %%
