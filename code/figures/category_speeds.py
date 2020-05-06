#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import seaborn as sns
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()


# Load the complex data set
data = pd.read_csv('../../data/compiled_estimate_categories.csv', comment='#')

# Define the division time window
div_range = [0.45, 0.55]

# Trim the data. 
data = data[(data['growth_rate_hr'] >= div_range[0]) 
        & (data['growth_rate_hr'] <= div_range[1]) & 
        (data['n_complex'] > 0)]

# Aggregate
summarized = data.groupby(['shorthand', 'name', 
                            'rate', 'rate_units'
                            ])['n_complex'].agg(
                            ('mean', 'min', 'max')).reset_index()

summarized['mean_rate_pool'] = summarized['mean'].values * summarized['rate'].values
summarized['min_rate_pool'] = summarized['min'].values * summarized['rate'].values
summarized['max_rate_pool'] = summarized['max'].values * summarized['rate'].values


#%%
# define the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(6, 3))

# Format the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_ylabel('time [sec]', fontsize=6)
ax.set_yscale('log')
ax.set_ylim(1, 1E4)

# Plot our stopwatch time. 
ax.hlines(6000, 0, 10, 'k', linestyle='--')

# Define the indices
indices = {'carbon':0,
           'atp': 1}

# Compute and populate categories. 

# Carbon Transport
N_sugars = 1E10 / 10
carbon_min, carbon_max, carbon_rate = summarized[
        summarized['shorthand']=='carbon_tport'][['min', 'max', 'rate']].values[0]
min_time = N_sugars / (carbon_max * carbon_rate)
max_time = N_sugars / (carbon_min * carbon_rate)
ax.vlines(indices['carbon'], min_time, max_time, lw=10, alpha=0.75)

# ATP synthesis
N_ATP = 4E9
atp_min, atp_max, atp_rate = summarized[
        summarized['shorthand']=='atp_synthase'][['min', 'max', 'rate']].values[0]
min_time = N_ATP / (atp_max * atp_rate)
max_time = N_ATP / (atp_min * atp_rate)
ax.vlines(indices['atp'], min_time, max_time, lw=10, alpha=0.75)

# RNA synthesis
sig_min, sig_max, sig_rate = summarized[
        summarized['shorthand']=='sigma70'][['min', 'max', 'rate']].values[0]
N_nts = 1E4 * 4500
min_time = N_rrnfi

# DNA synthesis
dna_min, dna_max, dna_rate = summarized[
        summarized['shorthand']=='atp_synthase'][['min', 'max', 'rate']].values[0]

# %%
#