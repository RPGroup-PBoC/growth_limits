#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import prot.viz 
import prot.estimate
colors = prot.viz.plotting_style()
constants = prot.estimate.load_constants()
dataset_colors = prot.viz.dataset_colors()

# Define some constants.
ribosomes = np.logspace(0, 7, 300)
r_aa = np.logspace(4, 8, 4)

def compute_elongation_rate(r_aa, R, Kd=1, V=1, t=1, Na=6.022E5, fa=1, rt_max=17.1):
    """
    Computes the physically meaningful root for the elongation rate. 
    """ 
    Kd *= Na
    root_term = np.sqrt(Kd**2 + 2* Kd * (rt_max * R * fa + r_aa) + (r_aa - rt_max * R * fa)**2)
    # a = -R * fa * t
    # b = r_aa * t + rt_max * R * fa * t + Kd * V * Na
    # c = -rt_max * r_aa * t
    # numer = -b + np.sqrt(b**2 - 4 * a *c)
    return (Kd + rt_max * R * fa + r_aa - root_term) / (2 * R * fa)
           
def compute_growth_rate(r_t, Naa, R, fa=1):
    """
    Given an elongation rate, compute the growth rate. 
    """
    return 3600 * r_t * R * fa * Naa**-1 

fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=200)
ax.set_xscale('log')
ax.set_ylim([0, 17.5])
ax.set_xlim([ribosomes[0], ribosomes[-1]])
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('ribosomes per µm$^{3}$', fontsize=6)
ax.set_ylabel('elongation rate [AA / s]', fontsize=6)

r_t = compute_elongation_rate(3E7, ribosomes)
ax.plot(ribosomes, r_t, '-', lw=1, color=colors['blue']) 

# plt.savefig('../../figures/r_aa_fixed_plots.svg')
# %%
