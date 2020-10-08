#%%
import numpy as np 
import pandas as pd 
import bokeh.io 
from bokeh.models import *
import bokeh.layouts
import prot.viz
colors, palette = prot.viz.bokeh_theme()
dataset_colors = prot.viz.dataset_colors()

# Load the data setc
cplx_desc = pd.read_csv('cplx_desc.csv')
data = pd.read_csv('../../data/')