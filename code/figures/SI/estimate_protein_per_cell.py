import numpy as np
import pandas as pd
from scipy import stats
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker
import prot.viz
import prot.size
colors, palette = prot.viz.bokeh_theme()
# dataset_colors = prot.viz.dataset_colors()
prot.viz.plotting_style()

# Exponential fit functions
from scipy.optimize import curve_fit
def func(x, a, c, d):
    return a*np.exp(-c*x)+d

##############################################
##############################################
# gather RNA/protein and DNA data from literature
##############################################
##############################################

# perform linear fit of RNA-ro-protein ratios
#  pairwise for < 0.7 hr-1 and >= 0.7 hr-1
dai_df = pd.read_csv('../../../data/dai2016_raw_data/dai2016_summary.csv')

dai_df_slow = dai_df[dai_df.growth_rate_hr < 0.7]
dai_df_fast = dai_df[dai_df.growth_rate_hr >= 0.7]

slope_dia_RP_A, intercept_dia_RP_A, r_value, p_value, std_err = stats.linregress(dai_df_slow.growth_rate_hr.values,
                                            dai_df_slow.RNA_P_ratio.values)
slope_dia_RP_B, intercept_dia_RP_B, r_value, p_value, std_err = stats.linregress(dai_df_fast.growth_rate_hr.values,
                                            dai_df_fast.RNA_P_ratio.values)

basan_df = pd.read_csv('../../../data/basan2015_raw_data/basan2015_data.csv')
popt_dna, pcov_dna = curve_fit(func, basan_df.growth_rate_hr.values, basan_df.dna_fg.values, p0=(1, 1e-6, 1))

##############################################
##############################################
# gather data from literature
##############################################
##############################################

# # Load the compiled data and grab schmidt growth rates
data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
data = data[data.dataset == 'schmidt_2016']

schmidt_gr = data.sort_values(by='growth_rate_hr', ascending = True).growth_rate_hr.unique() #np.array([lambda_dict[c] for c in data_condition])
# pred_vol_Si = (0.28 * np.exp(1.33  * schmidt_gr))
pred_vol_Si = prot.size.lambda2size(schmidt_gr)

# assume 1.1 g/ml mass density in cell, 30 % dry mass
# calculate total dry mass in fg
pred_drymass_Si = ((1.1*0.30)*(pred_vol_Si*1E-12)*1E15)

# calculate total DNA mass in fg
pred_dnamass_Si =  func(schmidt_gr, *popt_dna)

# Calculate total RNA + protein mass, assume DNA + RNA + protein == 90% dry mass
pred_RNA_protein_mass_Si = pred_drymass_Si*0.90 - pred_dnamass_Si

# predict RNA/Protein ratio:
pred_RNA_protein_ratio = np.append(\
                    (slope_dia_RP_A*schmidt_gr[:-2] + intercept_dia_RP_A),
                    (slope_dia_RP_B*schmidt_gr[-2:] + intercept_dia_RP_B))

# estimate total RNA mass per cell
pred_RNAmass_Si = (pred_RNA_protein_ratio * pred_RNA_protein_mass_Si) / (1+ pred_RNA_protein_ratio)

# estimate total protein mass per cell
pred_proteinmass_Si =  (pred_RNA_protein_mass_Si/ (1+pred_RNA_protein_ratio))

# ##############################################
##############################################
# Now plot!
##############################################
##############################################

fig, ax = plt.subplots(2, 2, figsize=(5, 4))
ax = ax.ravel()

# Plot DNA data from Basan et al.
# growth rates
x = np.linspace(0,2,100)

ax[0].plot(basan_df.growth_rate_hr.values, basan_df.dna_fg.values, 'o', ms=4, color=colors['light_purple'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'Basan et al 2015')
ax[0].plot(x, func(x, *popt_dna), '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[0].set_ylabel('DNA mass per cell [fg]', fontsize=6)
ax[0].legend(loc = 'upper left', fontsize=6)

# Plot RNA-to-protein data from Dai et al.
ax[1].plot(dai_df_slow.growth_rate_hr.values, dai_df_slow.RNA_P_ratio.values, 'o', ms=4, color=colors['pale_yellow'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'Dai et al 2016')
ax[1].plot(dai_df_fast.growth_rate_hr.values, dai_df_fast.RNA_P_ratio.values, 'o', ms=4, color=colors['pale_yellow'],
        markeredgewidth=0.5, markeredgecolor='k')

x = np.linspace(0,0.7,100)
ax[1].plot(x, (slope_dia_RP_A*x + intercept_dia_RP_A), '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)

x = np.linspace(0.7,2,100)
ax[1].plot(x, (slope_dia_RP_B*x + intercept_dia_RP_B), '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[1].set_ylabel('RNA: protein mass ratio', fontsize=6)
ax[1].legend(loc = 'upper left', fontsize=6)


# Plot predictions of DNA, RNA, and protein along with 'data points' for Schmidt
x = np.linspace(0,2,100)
vol = prot.size.lambda2size(schmidt_gr)
ax[2].plot(schmidt_gr, pred_dnamass_Si/vol, 'd', ms=4, color=colors['light_green'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'DNA')
ax[2].plot(schmidt_gr, pred_RNAmass_Si/vol, 's', ms=4, color=colors['light_green'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'RNA')
ax[2].plot(schmidt_gr, pred_proteinmass_Si/vol, 'o', ms=4, color=colors['light_green'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'protein')

ax[3].plot(schmidt_gr, pred_dnamass_Si, 'd', ms=4, color=colors['light_green'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'DNA')
ax[3].plot(schmidt_gr, pred_RNAmass_Si, 's', ms=4, color=colors['light_green'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'RNA')
ax[3].plot(schmidt_gr, pred_proteinmass_Si, 'o', ms=4, color=colors['light_green'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'protein')

# repeat calculations for smooth curve, predict RNA/Protein ratio <= 0.7 hr-1:
x = np.linspace(0,0.7,100)
vol = prot.size.lambda2size(x)
pred_drymass_Si = ((1.1*0.30)*(vol*1E-12)*1E15)
pred_dnamass_Si =  func(x, *popt_dna)
# Calculate total RNA + protein mass, assume DNA + RNA + protein == 90% dry mass
pred_RNA_protein_mass_Si = pred_drymass_Si*0.90 - pred_dnamass_Si
#  estimate R:P ratio
pred_RNA_protein_ratio = slope_dia_RP_A*x + intercept_dia_RP_A
# estimate total RNA mass per cell
pred_RNAmass_Si = (pred_RNA_protein_ratio * pred_RNA_protein_mass_Si) / (1+ pred_RNA_protein_ratio)
# estimate total protein mass per cell
pred_proteinmass_Si =  (pred_RNA_protein_mass_Si/ (1+pred_RNA_protein_ratio))

# concentration
ax[2].plot(x, func(x, *popt_dna)/vol, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[2].plot(x, pred_RNAmass_Si/vol, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[2].plot(x, pred_proteinmass_Si/vol, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
#  fg per cell
ax[3].plot(x, func(x, *popt_dna), '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[3].plot(x, pred_RNAmass_Si, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[3].plot(x, pred_proteinmass_Si, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)

# repeat calculations for smooth curve, predict RNA/Protein ratio >= 0.7 hr-1:
x = np.linspace(0.7, 2.0, 100)
vol = prot.size.lambda2size(x)
pred_drymass_Si = ((1.1*0.30)*(vol*1E-12)*1E15)
pred_dnamass_Si =  func(x, *popt_dna)
# Calculate total RNA + protein mass, assume DNA + RNA + protein == 90% dry mass
pred_RNA_protein_mass_Si = pred_drymass_Si*0.90 - pred_dnamass_Si
#  estimate R:P ratio
pred_RNA_protein_ratio = slope_dia_RP_B*x + intercept_dia_RP_B
# estimate total RNA mass per cell
pred_RNAmass_Si = (pred_RNA_protein_ratio * pred_RNA_protein_mass_Si) / (1+ pred_RNA_protein_ratio)
# estimate total protein mass per cell
pred_proteinmass_Si =  (pred_RNA_protein_mass_Si/ (1+pred_RNA_protein_ratio))

# concentration
ax[2].plot(x, func(x, *popt_dna)/vol, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[2].plot(x, pred_RNAmass_Si/vol, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[2].plot(x, pred_proteinmass_Si/vol, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
#  fg per cell
ax[3].plot(x, func(x, *popt_dna), '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[3].plot(x, pred_RNAmass_Si, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)
ax[3].plot(x, pred_proteinmass_Si, '-', alpha = 0.6,
            lw = 0.5, color = 'k', zorder = 0)


ax[2].set_ylabel('concentration [fg/fL]', fontsize=6)
ax[2].set_ylim(0,270)

ax[3].set_ylabel('protein mass per cell [fg]', fontsize=6)
ax[3].legend(loc = 'upper left', fontsize=6)
ax[3].set_ylim(0,850)

# for ax_ in ax:
#     for c, d in size_df.groupby(['dataset_name']):
#         ax_.plot(d['growth_rate_hr'].values, d['volume_um3'].values, 'o', ms=4, color=colordict[c],
#                 markeredgewidth=0.5, markeredgecolor='k', label = c)
#         ax_.set_ylabel(r'average cell size [$\mu$m$^3$]', fontsize=6)
#
#         ax_.plot(x, 0.28 * np.exp(1.33  * x), '-', alpha = 0.6,
#                     lw = 0.5, color = 'k', zorder = 0)
#         ax_.plot(x, 0.27 * 2**(1.1  * x / np.log(2)), '-', alpha = 0.6,
#                     lw = 0.5, color = 'k', zorder = 0)
#         ax_.plot(x[:-30], func(x[:-30], *popt_vol), '-', alpha = 0.6,
#                     lw = 0.5, color = 'k', zorder = 0)

for ax_ in ax:
    ax_.xaxis.set_tick_params(labelsize=5)
    ax_.yaxis.set_tick_params(labelsize=5)
    ax_.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    ax_.set_xlim([0, 2.0])
#
# ax[1].set_yscale('log')
# ax[1].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# ax[1].get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
#
# ax[1].set_yticks([0.3, 0.5, 1, 2, 3, 5])
#
# handles, labels = ax[0].get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# ax[0].legend(by_label.values(), by_label.keys(), loc = 'upper left', fontsize = 6)
#

plt.tight_layout()
fig.savefig('../../../figures/predict_protein_per_cell.pdf', bbox_inches='tight')
