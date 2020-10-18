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
div_range = [0.48, 0.52]

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

# Define a function that will report the min and max bounds
def bounds(process, N, df=summarized):
        min_num, max_num, rate = summarized[summarized['shorthand']==process][
                                   ['min', 'max', 'rate']].values[0]
        min_time = N / (max_num * rate)
        max_time = N / (min_num * rate)
        return [min_time, max_time]

#%%
# define the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 4))

# Format the axes
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('time [sec]', fontsize=6)
ax.set_xscale('log')
ax.set_xlim(90, 7E3)
ax.set_xticks([100, 250, 1000, 2500])
ax.set_xticklabels(['10$^{2}$', '2.5$\times$ 10$^2$', '10$^3$', '2.5$\times$10$^3$'])
ax.set_ylim([-0.5,6.5])
# Plot our stopwatch time. 
ax.vlines(6000, -1, 7, 'k', linestyle='--')

# Define the indices
viridis = sns.color_palette('viridis', n_colors=10)
indices = {'protein': 0,
           'trna':1,
            'atp': 2,
           'dntp': 3,
           'sigma70': 4,
           'carbon':5,
           'dnap': 6}
labels = ['protein synthesis', 'tRNA charging', 'ATP syntheis',
          'dNTP synthesis', 'rRNA synthesis', 'carbohydrate transport',
          'DNA synthesis']
ax.set_yticks(list(indices.values()))
ax.set_yticklabels(labels)

process_colors = {'carbon':viridis[indices['carbon']],
                  'atp': viridis[indices['atp']],
                  'sigma70': viridis[indices['sigma70']],
                  'dnap': viridis[indices['dnap']],
                  'dntp': viridis[indices['dntp']],
                  'protein': viridis[indices['protein']],
                  'trna':viridis[indices['trna']]}


# Compute and populate categories. 

# Carbon Transport
N_sugars = 1E10 / 10
min_time, max_time = bounds('carbon_tport', N_sugars)
ax.hlines(indices['carbon'], min_time, max_time, lw=2, alpha=0.75,
          color=process_colors['carbon'])
ax.plot(min_time, indices['carbon'], 'o', ms=4.5, color=process_colors['carbon'])
ax.plot(max_time, indices['carbon'], 'o', ms=4.5, color=process_colors['carbon'])

# ATP synthesis
N_ATP = 4E9
min_time, max_time = bounds('atp_synthase', N_ATP)
ax.hlines(indices['atp'], min_time, max_time, lw=2, alpha=0.75, 
          color=process_colors['atp'])
ax.plot(min_time, indices['atp'], 'o', ms=4.5, color=process_colors['atp'])
ax.plot(max_time, indices['atp'], 'o', ms=4.5, color=process_colors['atp'])


# RNA synthesis
N_nts = 1E4 * 4500
min_time, max_time = bounds('sigma70', N_nts)
ax.hlines(indices['sigma70'], min_time, max_time, lw=2, alpha=0.75,
          color=process_colors['sigma70'])
ax.plot(min_time, indices['sigma70'], 'o', ms=4.5, color=process_colors['sigma70'])
ax.plot(max_time, indices['sigma70'], 'o', ms=4.5, color=process_colors['sigma70'])


# DNA synthesis
N_bp = 5E6
min_time, max_time = bounds('dnap', N_bp)
ax.hlines(indices['dnap'], min_time, max_time, lw=2, alpha=0.75,
        color=process_colors['dnap'])
ax.plot(min_time, indices['dnap'], 'o', ms=4.5, color=process_colors['dnap'])
ax.plot(max_time, indices['dnap'], 'o', ms=4.5, color=process_colors['dnap'])


N_dntp = 5E6
min_time, max_time = bounds('dntp', N_dntp)
ax.hlines(indices['dntp'], min_time, max_time, lw=2, alpha=0.75,
        color=process_colors['dntp'])
ax.plot(min_time, indices['dntp'], 'o', ms=4.5, color=process_colors['dntp'])
ax.plot(max_time, indices['dntp'], 'o', ms=4.5, color=process_colors['dntp'])


# Protein synthesis
N_AA = 1E9
min_time, max_time = bounds('ribosome', N_AA)
ax.hlines(indices['protein'], min_time, max_time, lw=2, alpha=0.75,
        color=process_colors['protein'])
ax.plot(min_time, indices['protein'], 'o', ms=4.5, color=process_colors['protein'])
ax.plot(max_time, indices['protein'], 'o', ms=4.5, color=process_colors['protein'])

N_tRNA = 1E9
min_time, max_time = bounds('trna', N_tRNA)
ax.hlines(indices['trna'], min_time, max_time, lw=2, alpha=0.75,
        color=process_colors['trna'])
ax.plot(min_time, indices['trna'], 'o', ms=4.5, color=process_colors['trna'])
ax.plot(max_time, indices['trna'], 'o', ms=4.5, color=process_colors['trna'])
plt.savefig('../../figures/synthesis_times.svg', bbox_inches='tight')




# %%
#