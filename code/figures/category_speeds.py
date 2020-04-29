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
div_range = [0.4, 0.6]

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


# Define the numbers needed
num_needed = {'rnap': 2E3 * 1E3,
              'dnap': 5E6,
              'glucose_tport':  1E10,
              'atp_synthase': 1E10,
              'ribosome': 2E4 * 1E4,
              'fas': 5E7,
              'ndhI': 3E9,
              'trna': 5E6}

for k, v in num_needed.items():
    summarized.loc[summarized['shorthand']==k, 'num_needed'] = v

# Compute the time it would take to  achieve estimate
summarized['mean_time'] = summarized['num_needed'].values / summarized['mean_rate_pool'].values
summarized['min_time'] = summarized['num_needed'].values / summarized['max_rate_pool'].values
summarized['max_time'] = summarized['num_needed'].values / summarized['min_rate_pool'].values
summarized.sort_values('mean_time', inplace=True)
summarized['location'] = np.arange(0, len(summarized), 1)
summarized['color'] = sns.color_palette('viridis', n_colors=len(summarized) + 3)[:len(summarized)]
# %%
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(4, 6))
ax.set_yticks(summarized['location'].values)
ax.set_yticklabels(summarized['name'].values)
ax.set_xlabel('time [s]')
ax.set_ylim([-1, 8])
ax.vlines(6E3, -1, 9, 'k', linestyle='dashed')
ax.set_xscale('log')

for g, d in summarized.groupby(['name', 'location', 'color']):
    ax.plot(d['min_time'], g[1], 'o', color=g[2])
    ax.plot(d['max_time'], g[1], 'o', color=g[2])
    ax.plot([d['min_time'].values[0], d['max_time'].values[0]], [g[1], g[1]], '-', color=g[2])
plt.savefig('../figures/observed_time_categories.pdf', bbox_inches='tight')
# %%


# %%
