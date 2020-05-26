#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Load and restrict the data
data = pd.read_csv('../../data/compiled_estimate_categories.csv')
ribo = data[data['shorthand']=='ribosome']

# Define the color palettes. 
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}


# Define some parameters
growth_rate = np.linspace(0, 2, 200)
t_double = growth_rate**-1 * 60 * 60
DENSITY = 1.1 # in pg per fL
PROT_FRAC = 0.15 # Fraction of dry mass that is protein
AVG_PROT_MASS = 30E3 / (6E23 * 1E-12) # average protein mass in pg
AVG_PROT_LENGTH = 300 # in Amino Acids
AVG_AA_MASS = 110 / (6E23 * 1E-12)
TRANSLATION_RATE = 17.1 # in AA per sec
volume =  0.28 * np.exp(1.33  * growth_rate) # From Si et al. 2017
a = -0.89
c = 4.60
d = 0.922
active_frac = a*np.exp(-c*growth_rate)+d


# Compute the predicted scaling. 
prediction = (DENSITY * volume * PROT_FRAC) /\
             (AVG_AA_MASS * active_frac * TRANSLATION_RATE * t_double)



# %%
fig, ax = plt.subplots(1, 1, figsize=(6,4))
# Format the axes
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribosomes per cell', fontsize=6)

# Plot the prediction
ax.plot(growth_rate, prediction, 'k--', lw=0.75, label='predicted trend')


# Plot the measurements
for g, d in ribo.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, 
            color=dataset_colors[g[0]], markeredgecolor='k', markeredgewidth=0.5,
            alpha=0.75, label=g[1])

ax.legend(fontsize=6)
plt.savefig('../../figures/ribosome_trend.pdf')
# %%
