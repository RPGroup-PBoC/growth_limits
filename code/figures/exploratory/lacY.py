#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import prot.viz
# prot.viz.plotting_style()
data = pd.read_csv('../../../data/compiled_annotated_complexes.csv')
data = data[data['gene_name'].str.contains('lac')]

data.head()

# %%
