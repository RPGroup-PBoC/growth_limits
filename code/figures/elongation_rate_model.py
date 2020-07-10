#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import prot.viz 
import prot.estimate
colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

# Load the estimate data 
data = pd.read_csv('../../data/compiled_estimate_categories.csv')


# Compute concentrations per µm^3
data['conc'] = data['n_complex'].values / data['volume']

# COmpute the averate tRNA synthase concentration
tRNA = data[data['shorthand']=='trna']['conc']
mean_tRNA = np.mean(tRNA.values)
k_trna = 20 # in charging events 

# As a function of growth rate, compute the amino acid supply. 
r_aa =  mean_tRNA * k_trna

# Compute the ribosome range
ribosomes = np.logspace(3, 5, 300)
kd = 5
avo = 6.022E5
N_AA = 1.1 * 0.15 / (110 / 6E11) # Number of amino acids per µm^3 (approx)
N_AA_ribosomes = N_AA + (7500 * ribosomes) 

# Compute the elongation rate. 
theta = (kd * avo) / (np.log(2) * N_AA_ribosomes)
rt_max = 17.1

# Define a function to compute the predicted growth rate
sqrt_term = np.sqrt(4 * theta * rt_max * r_aa * ribosomes + (rt_max * ribosomes)**2 - 2 * rt_max * r_aa * ribosomes + r_aa**2)
r_t = (rt_max * ribosomes + r_aa - sqrt_term) / (2 * ribosomes * (1 - theta))

predicted_growth_rate = 3600 * (r_t * ribosomes) / N_AA_ribosomes

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.set_xscale('log')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('ribosomes per µm$^3$', fontsize=6)
ax.set_ylabel('growth rate λ [hr$^{-1}$]', fontsize=6)
# ax.set_ylim([0 , 2])
# ax.set_xlim([1E3, 1E5])

for g, d in data[data['shorthand']=='ribosome'].groupby(['dataset', 'dataset_name']):
    ax.plot(d['conc'], d['growth_rate_hr'], 'o', ms=5,
            color=dataset_colors[g[0]], label=g[1], markeredgewidth=0.5,
            markeredgecolor='k')

ax.plot(ribosomes, predicted_growth_rate, 'k-')
# %%fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.set_yscale('log')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('tRNA synthase per µm$^3$', fontsize=6)
# ax.set_xlim([0 , 2])
# ax.set_ylim([1E3, 1E5])

for g, d in data[data['shorthand']=='trna'].groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'].values , d['n_complex'].values / d['volume'], 'o', ms=5,
            color=dataset_colors[g[0]], label=g[1], markeredgewidth=0.5,
            markeredgecolor='k')



# %%
