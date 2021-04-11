#%%
import numpy as np
import pandas as pd
import prot.viz
import matplotlib.pyplot as plt
colors = prot.viz.plotting_style()

def elongation_rate(R, r_AA, V=1E-15, K_D=5E-3, t=1, f_a = 1, rt_max=17.1):

    """
    Computes the average elongation rate of actively translating ribosomes
    given the model parameters.

    Parameters
    ----------
    R : array_like, float or int
        The number of ribosomes in the desired cell volume. Must be greater than 0
    r_AA: array_like, float or int
        The amino acid supply rate to the volume in units of AA per ribosome
    V : array_like, float or int, optional
        The unit volume to consider in units of L. Default value is 1 fL.
    K_D: array_like, float, optional
        The effective dissociation constant of aa-tRNAs to the ribosome in units
        of M. Default value is 5 mM
    t : float, optional
        The timescale in units of seconds. Default value is 1 second.
    f_a: float, optional
        The fraction of the ribosome pool that is actively translating. Default value
        is 1, implying all ribosomes are translating.
    rt_max: float, optional
        The maximal elongation rate in units of AA per second per ribosome. Default
        value is 17.1 AA / sec / ribosome, which is the maximal elongation rate
        in E. coli

    Returns
    -------
    r_t: array_like
        The average elongation rate of all actively elongating ribosomes.
    """

    N_A = 6.022E23 # Avogadro's number
    prefix = K_D * N_A * V + R * f_a * rt_max * t + r_AA * t
    sqrt_a = (K_D * N_A * V)**2 + (r_AA * t)**2 + (R * f_a * rt_max * t)**2
    sqrt_b = 2 * (K_D * N_A * R * V * f_a * rt_max* t + K_D * N_A * V * r_AA * t -\
                 R * f_a * rt_max * r_AA * t**2)
    denom = 2 * R * f_a * t
    r_t = (prefix - np.sqrt(sqrt_a + sqrt_b)) / denom
    return r_t

def growth_rate_indep(phi_R, V, r_aa, f_a=1):
        L_R = 7459
        N_pep = (V * 1.1E-12 * 0.3 * 0.5 * 6.022E23 / 110)
        R = phi_R * N_pep / L_R
        R_conc = R/V
        r_aa_conc = r_aa/V
        r_t = elongation_rate(R_conc, r_aa_conc)
        return 3600 * r_t * R * f_a / N_pep


phi_R_range = np.linspace(0.05, 0.5, 200)
V_range = np.linspace(0.3, 2, 200)
rAA_range = [1E4, 5E5, 5E6, 1E8]
phi_R, V = np.meshgrid(phi_R_range, V_range)
heat_maps = []

for i, r in enumerate(rAA_range):
    out = growth_rate_indep(phi_R, V, r)
    heat_maps.append(out)

fig, ax = plt.subplots(2, 2, figsize=(6, 6))
for i, a in enumerate(ax.ravel()):
    pcm = a.imshow(heat_maps[i], origin='lower', vmin=0, vmax=3, cmap='viridis')
    conts = a.contour(np.arange(0, 200), np.arange(0, 200),
                        heat_maps[i], zorder=1000, colors='white')
    # clabel = a.clabel(conts)
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlabel('$\Phi_R$', fontsize=8)
    a.set_ylabel('volume [µm$^3$]', fontsize=8)
    a.set_title('$r_{AA}$ = 10$^{' + f'{int(np.log10(rAA_range[i]))}' + '}$ AA / sec / cell',
                fontsize=8, backgroundcolor=colors['pale_yellow'], y=1.08)
    a.set_yticks([0, 50, 100, 150, 200])
    a.set_xticks([0, 50, 100, 150, 200])
    a.set_yticklabels([str(np.round(V_range[i], decimals=1)) for i in [0, 49, 99, 149, 199]])
    a.set_xticklabels([str(np.round(phi_R_range[i], decimals=1)) for i in [0, 49, 99, 149, 199]])
fig.subplots_adjust(hspace=0.5, wspace=0.01)
bar = fig.colorbar(pcm, ax=ax, orientation='horizontal',
                    fraction=0.1, shrink=0.8, pad=0.08)
fig.text(0.12, 0.915, '(A)', fontsize=8)
fig.text(0.5, 0.915, '(B)', fontsize=8)
fig.text(0.12, 0.545, '(C)', fontsize=8)
fig.text(0.5, 0.545, '(D)', fontsize=8)
bar.ax.tick_params(labelsize=6)
bar.set_label('growth rate [hr$^{-1}$]', fontsize=8)
plt.savefig('../../figures/figS12_RV_heatmaps.pdf', bbox_inches='tight')
# bar.set_tickparams(labelsize=6)
# bar.set_title('growth rate [hr$^{-1}$]', fontsize=8)
# plt.tight_layout()

# %%
