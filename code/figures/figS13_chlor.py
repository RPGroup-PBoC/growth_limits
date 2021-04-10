#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size as size
prot.viz.plotting_style()
colors, palette = prot.viz.bokeh_theme()
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def elongation_rate(R, r_AA, V, K_D, t, f_a, rt_max):

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
        The unit volume to consider in units of L.
    K_D: array_like, float, optional
        The effective dissociation constant of aa-tRNAs to the ribosome in units
        of M.
    t : float, optional
        The timescale in units of seconds.
    f_a: float, optional
        The fraction of the ribosome pool that is actively translating. Default value
        is 1, implying all ribosomes are translating.
    rt_max: float, optional
        The maximal elongation rate in units of AA per second per ribosome.

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


def compute_growth_rate(r_t, Naa, R, fa=1):
    """
    Given an elongation rate, compute the growth rate.
    """
    return 3600 * r_t * R * fa * Naa**-1


# load in the data from Dai et al. 2016:
dai_Cm_df = pd.read_csv('../../data/dai2016_raw_data/dai2016_chlor_data_summary.csv')
dai_nut_df = pd.read_csv('../../data/dai2016_raw_data/dai2016_nutrient_data_summary.csv')

color_dict = dict(zip(np.append(dai_Cm_df.condition.unique(), dai_nut_df.condition.unique()), palette))


# color map for varying fraction of active ribosomes:
cm = plt.cm.get_cmap('viridis')

# initialize plot
fig, ax = plt.subplots(1, 2, figsize = (8,4))

# we need to 'infer' best values for r_aa for each condition without
# antibiotic; we'll assume this doesn't change with the addition of chlor.
# place values in  df_raa_est
df_raa_est = pd.DataFrame()

for gr, d in dai_nut_df.groupby(['growth_rate_hr', 'condition']):
    r_aa = np.logspace(-1,10, 1000)
    # rt = compute_elongation_rate(r_aa, d['R'].values[0], Kd=1E-1, V=d['V'].values[0]*1E-15, #Kd=1E-5
    #                              t=60*60*(np.log(2)/d['growth_rate_hr'].values[0]), #t=d['growth_rate_hr'].values[0],
    #                     Na=6.022E23, fa=d['f_a'].values[0], rt_max=18)

    rt  = elongation_rate( d['R'].values[0], r_aa, V=d['V'].values[0]*1E-15, K_D=1E-1,
            t=60*60*(np.log(2)/d['growth_rate_hr'].values[0]), f_a = d['f_a'].values[0], rt_max=18)


    rt_compare = list(np.abs(rt - d['Translational elongation rate (aa/s)'].values[0]))

    minpos = rt_compare.index(min(rt_compare))


    df_raa_est_list = {'condition' :  gr[1],
                      'growth_rate_hr' : gr[0],
                      'raa_bestfit' : r_aa[minpos],
                      'rt' : d['Translational elongation rate (aa/s)'].values[0],
                      'f_a' : d['f_a'].values[0],
                      'rt_fit' : rt[minpos],
                      'V' : d['V'].values[0],
                      'R' : d['R'].values[0],
                      'Naa' : d['Naa'].values[0],
                      'Phi_R' : d['Phi_R'].values[0]}
    df_raa_est = df_raa_est.append(df_raa_est_list,
                                  ignore_index = True)

    ax[0].scatter(gr[0], rt[minpos], zorder=10,
               color = 'k', alpha = 0.5)

    # also plot effect as fcn of f_a
    f_a = np.linspace(0,d['f_a'].values[0],100)
    # rt_fa = compute_elongation_rate(r_aa[minpos], d['R'].values[0], Kd=1E-1, V=d['V'].values[0]*1E-15, #Kd=1E-5
    #                          t=60*60*(np.log(2)/d['growth_rate_hr'].values[0]), #t=d['growth_rate_hr'].values[0],
    #                 Na=6.022E23, fa=f_a, rt_max=18)

    rt_fa  = elongation_rate( d['R'].values[0], r_aa[minpos], V=d['V'].values[0]*1E-15, K_D=1E-1,
            t=60*60*(np.log(2)/d['growth_rate_hr'].values[0]), f_a = f_a, rt_max=18)


    gr_fa = compute_growth_rate(rt_fa, d['Naa'].values[0], d['R'].values[0], f_a)

    sc = ax[0].scatter(gr_fa, rt_fa, c=f_a, cmap=cm, vmin = 0,  vmax =  1, s = 4, alpha = 0.5)

ax[0].set_ylim(8,17.5)
ax[0].set_xlim(0,2)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize = 14)
ax[0].set_ylabel('elongation rate [aa/s]', fontsize = 14)
cbaxes = inset_axes(ax[0], width="30%", height="4%", loc=4)
plt.colorbar(sc, cax=cbaxes, ticks=[0.,1], orientation='horizontal')

cbaxes.set_xlabel('active ribosomal fraction $f_a$')



for cond, d in dai_Cm_df.groupby('condition'):
    df_raa_est_ = df_raa_est[df_raa_est.condition == cond]

    rt = 1/((1/20) + 1/(6.4*31*d['RNA_P_ratio']))
    gr = compute_growth_rate(rt, df_raa_est_.Naa.values[0], df_raa_est_.R.values[0], fa=d['f_a'].values)

    ax[1].plot(gr, rt, zorder=10, ls = '--',
               color =  color_dict[cond])

    ax[1].scatter(gr, d['Translational elongation rate (aa/s)'], zorder = 10,
              color =  color_dict[cond])


ax[1].set_ylim(8,17.5)
ax[1].set_xlim(0,2)
ax[1].set_xlabel('growth rate', fontsize = 14)
ax[1].set_ylabel('elongation rate [aa/s]', fontsize = 14)

plt.tight_layout()

plt.savefig('../../figures/figS13_Cm_data.pdf')
