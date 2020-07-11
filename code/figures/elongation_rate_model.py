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
vol = 1 # in cubic microns
ribosomes = np.logspace(3, 6, 300)  / vol # Ribosomes per volume
r_aa = np.logspace(4, 7, 4) / vol # Amino acid supply rate per second per volume
Na = 6.022E5 # For conversion of concentrations to mM 
Kd = 5 # Dissociation constant in mM for abstracted AA to ribosome. 
Lr = 7500 # Length of a ribosome in amino acids. 
m_aa = 110 / 6E11 # Mass of amino acid in picograms
rho = 1.1 # Cell density in pg per fL
m_dry = vol * 1.1 * 0.3 # dry mass in pg
theta_prot = 0.5 # Fraction of dry mass that is protein
background_N_aa =  theta_prot * m_dry / m_aa # Background nummber of amino acids
N_aa = background_N_aa + ribosomes * Lr # Increasing number of amino acis based on addition of ribosomes.

# Define a funciton to compute the elongation rate. 
def compute_elongation_rate(Kd, r_aa, R, V=vol, N_aa=N_aa, Na=Na, fa=1, rt_max=17.1):
    """
    Computes the physically meaningful root for the elongation rate. 
    """
    theta = (Kd * V * Na) / (np.log(2) * N_aa)
    root_term = np.sqrt(4 * theta * rt_max * r_aa * R * fa +\
                         (rt_max * R * fa)**2 - 2 * rt_max * r_aa * R * fa +\
                          r_aa**2)
    r_t = (rt_max * R * fa + r_aa - root_term) / (2 * R * fa * (1 - theta))
    return r_t


def compute_growth_rate(r_t, Naa, R, fa=1):
    """
    Given an elongation rate, compute the growth rate. 
    """
    return 3600 * r_t * R * fa * Naa**-1 


elongation_colors = sns.color_palette('viridis', n_colors=len(r_aa) + 3) 
growth_rate_colors = sns.color_palette('viridis', n_colors=len(r_aa) + 3) 

# set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
for a in ax:
    a.set_xscale('log')
    a.set_xlim([ribosomes[0], ribosomes[-1]])
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlabel('ribosomes per µm$^{3}$', fontsize=6)

ax[0].set_ylabel('elongation rate [AA / s]', fontsize=6)
ax[1].set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)

for i, r in enumerate(r_aa):
    r_t = compute_elongation_rate(Kd, r, ribosomes)
    gr = compute_growth_rate(r_t, N_aa, ribosomes)
    ax[0].plot(ribosomes, r_t, '-', lw=1, color=elongation_colors[i], 
                label=r'10$^{%s}$ ' %int(np.log10(r)))
    ax[1].plot(ribosomes, gr, '-', lw=1, color=growth_rate_colors[i], 
                label=r'$r_{AA}$ = 10$^{%s}$' %int(np.log10(r)))

leg = ax[1].legend(fontsize=6, title='r$_AA$ [AA/s•µm$^3$]')
leg.get_title().set_fontsize(6)
plt.savefig('../../figures/r_aa_variation_plots.svg')
# %%
# set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
for a in ax:
    a.set_xscale('log')
    a.set_xlim([ribosomes[0], ribosomes[-1]])
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlabel('ribosomes per µm$^{3}$', fontsize=6)

ax[0].set_ylabel('elongation rate [AA / s]', fontsize=6)
ax[1].set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)

r_t = compute_elongation_rate(Kd, r_aa[2], ribosomes)
gr = compute_growth_rate(r_t, N_aa, ribosomes)
ax[0].plot(ribosomes, r_t, '-', lw=1, color=elongation_colors[2], 
            label=r'10$^{%s}$ ' %int(np.log10(r)))
ax[1].plot(ribosomes, gr, '-', lw=1, color=growth_rate_colors[2], 
            label=r'$r_{AA}$ = 10$^{%s}$' %int(np.log10(r)))

leg = ax[1].legend(fontsize=6, title='r$_AA$ [AA/s•µm$^3$]')
leg.get_title().set_fontsize(6)
plt.savefig('../../figures/r_aa_fixed_plots.svg')


# %%
