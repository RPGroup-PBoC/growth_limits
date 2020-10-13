#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import prot.viz
import prot.estimate
from mpl_toolkits.mplot3d import Axes3D

colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

# Define some constants.
ribosomes = np.logspace(1, 7, 300)


def compute_elongation_rate(r_aa, R, Kd=5E-3, V=1E-15, t=1, Na=6.022E23, fa=1, rt_max=17.1):
    """
    Computes the physically meaningful root for the elongation rate.
    """
    a = -R * fa * t
    b = r_aa * t + rt_max * R * fa * t + Kd * V * Na
    c = -rt_max * r_aa * t
    numer = -b + np.sqrt(b**2 - 4 * a *c)
    return numer / (2 * a)

def compute_growth_rate(r_t, Naa, R, fa=1):
    """
    Given an elongation rate, compute the growth rate.
    """
    return 3600 * r_t * R * fa * Naa**-1

fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), dpi=200)
ax.set_xscale('log')
# ax.set_ylim([0, 17.5])
ax.set_xlim([ribosomes[0], ribosomes[-1]])
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('ribosomes per µm$^{3}$', fontsize=6)
ax.set_ylabel('elongation rate [AA / s]', fontsize=6)

r_aa = 5E6
r_t = compute_elongation_rate(r_aa, ribosomes)
ax.plot(ribosomes[ribosomes < 3.3E5], r_t[ribosomes < 3.3E5], '-', lw=1, color=colors['blue'])
ax.plot(ribosomes[ribosomes > 3.3E5], r_t[ribosomes > 3.3E5], '--', lw=1, color=colors['blue'])

# Add points for illustration

ax.plot(500, compute_elongation_rate(r_aa, 500), 'o', ms=7, color=colors['light_blue'], markeredgewidth=0.25,
        markeredgecolor='k')
ax.plot(1E5, compute_elongation_rate(r_aa, 1E5), 'o', ms=7, color=colors['light_blue'], markeredgewidth=0.25,
        markeredgecolor='k')
ax.plot(1E6, compute_elongation_rate(r_aa, 1E6), 'o', ms=7, color=colors['light_blue'], markeredgewidth=0.25,
        markeredgecolor='k')
plt.savefig('../../figures/fig12_rt_vs_R.svg')


# %%
# Plot the growth rate vs ribosome copy number per cell for different r_aa
ribosomes = np.logspace(3.5, 5.7, 200)
# Compute the total number of amino acids as a function of ribosomes from the fit
NAA = np.exp(16.1) * ribosomes**0.5
vol = np.exp(-5.8) * ribosomes**0.6

# Compute the elongation rate
r_aa = np.logspace(4, 8.3, 10)# [1E5, 5E5, 1E6, 5E6, 1E7, 3E7]


fig, ax = plt.subplots(1, 1, figsize=(3, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('ribosomes per cell', fontsize=6)
ax.set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_xscale('log')
ax.set_xlim([ribosomes[0], ribosomes[-1]])
# ax.set_yscale('log')


# Set the color palette
viridis = sns.color_palette('viridis', n_colors=len(r_aa) + 1)
for i, r in enumerate(r_aa):
    r_t = compute_elongation_rate(r/vol, ribosomes, V=vol * 1E-15)
    growth = compute_growth_rate(r_t, NAA, ribosomes)
    ax.plot(ribosomes, growth, '-', lw=1, color=viridis[i],
            label='10$^{%s}$' %np.round(np.log10(r), decimals=1))
    ind = np.argmax(growth)
    # ax.plot(ribosomes[ind], growth[ind], 'v', ms=5, color=viridis[i],
        #    markeredgewidth=0.25, markeredgecolor='k', label='__nolegend__')

# ax.plot([], [], 'v', ms=5, color='grey', markeredgecolor='k',
    #    markeredgewidth=0.25, label='maximum growth rate')
leg = ax.legend(title=r'r$_{AA}$ [$\frac{AA}{s\cdot \mu m^3}$]',
               fontsize=6, bbox_to_anchor=(1,1))
leg.get_title().set_fontsize(6)
plt.savefig('../../figures/fig12_lambda_v_R.svg', bbox_inches='tight')


#%% Connection with data
data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# Compute the concentration of the tRNA synthases per µm^3 and use this as a
# measure of r_aa

data['conc'] = data['n_complex'].values / data['volume']
trna = data[data['shorthand']=='trna']['conc'].agg(('mean', 'std'))
r_aa_mean  = trna['mean'] * 00
# r_aa_max = (trna['mean'] + trna['std']) * 500
# r_aa_min = (trna['mean'] - trna['std']) * 500

R = data[data['shorthand']=='ribosome']
# Compute the elongation rate
ribo_range = np.logspace(3, 5.5, 200)
# Compute the total number of amino acids as a function of ribosomes from the fit
NAA_range = np.exp(16.1) * ribo_range**0.5
vol_range = np.exp(-5.8) * ribo_range**0.6
# rt_mean_data = compute_elongation_rate(3E6, R['conc'])
rt_mean_range = compute_elongation_rate(, ribo_range/vol_range, fa=0.5)
# rt_max_data = compute_elongation_rate(r_aa_max, R['conc'])
# rt_max_range = compute_elongation_rate(r_aa_max, ribo_range/vol_range)
# rt_min_data = compute_elongation_rate(r_aa_min, R['conc'])
# rt_min_range = compute_elongation_rate(r_aa_min, ribo_range/vol_range)

# compute the growth rate.
gr_mean_range = compute_growth_rate(rt_mean_range, NAA_range, ribo_range, fa=0.5)

fig, ax = plt.subplots(1,1)
ax.set_xscale('log')
ax.plot(ribo_range, gr_mean_range, 'k-')
ax.plot(R['n_complex'], R['growth_rate_hr'], 'o')

#%%
R['conc']

#%%
ribosomes = np.logspace(3, 4, 200)
fig, ax = plt.subplots(3, 1, figsize=(5, 6), sharex=True, dpi=100)
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xscale('log')
    # a.set_yscale('log')
ax[1].set_xlabel('ribosomes per cell', fontsize=6)
ax[0].set_ylabel('elongation rate [AA / s]', fontsize=6)
ax[1].set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)

# Set the limits
# ax[0].set_ylim([10, 12])
# ax[1].set_ylim([1E6, 1E7])

rt_low = compute_elongation_rate(3E3, ribosomes/vol, V=vol* 1E-15)
gr_high = compute_growth_rate(rt_high, NAA, ribosomes)
gr_low = compute_growth_rate(rt_low, NAA, ribosomes)
rt_reg = compute_elongation_rate(3E3, ribosomes/vol ,fa=0.1, V=vol * 1E-15)
gr_reg = compute_growth_rate(rt_reg, NAA, ribosomes, fa=0.1)
pref = (6.022E5 * vol)
# ax[0].plot(ribosomes, rt_high, 'k-', lw=0.75, label=r'r$_{AA}\approx 5 = 10^6 \frac{AA}{s \times cell}$ ; $f_a = 1$')
# ax[1].plot(ribosomes, gr_high, 'k-', lw=0.75, label='__nolegend__')
ax[0].plot(ribosomes, rt_low, '-', color=colors['red'], lw=0.75, label=r'r$_{AA}= 3 \times 10^6 \frac{AA}{s \times cell}$ ; $f_a = 1$')
ax[1].plot(ribosomes, gr_low, '-', color=colors['red'], lw=0.75, label='__nolegend__')
ax[0].plot(ribosomes, rt_reg, '-', color=colors['blue'], lw=0.75, label=r'r$_{AA} = 3 \times 10^6 \frac{AA}{s \times cell}$ ; $f_a = 0.5$')
# ax[1].plot(ribosomes, 1E-6 * (5E6 - (r_t * ribosomes)), 'k-', lw=0.75)
ax[2].plot(ribosomes, pref**-1 * (3E3 - (rt_low * ribosomes)), '-', lw=0.75, color=colors['red'])
ax[2].plot(ribosomes, pref**-1 * (3E3 - (rt_reg * ribosomes * 0.1)), '-', lw=0.75, color=colors['blue'])
ax[1].plot(ribosomes, gr_reg, '-', color=colors['blue'],lw=0.75, label='__nolegend__')
ax[0].legend(fontsize=6, bbox_to_anchor=(1,1))
# plt.savefig('../../figures/fa_v_R.pdf')
# %%


# %%
