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

# Exponential fit function
from scipy.optimize import curve_fit
def func(x, a, c, d):
    return a*np.exp(-c*x)+d

basan_df = pd.read_csv('../../../data/basan2015_raw_data/basan2015_data.csv')


##############################################
##############################################
# Perform exponential fit of Basan
##############################################
##############################################

popt_fg, pcov_fg = curve_fit(func, basan_df.growth_rate_hr.values, basan_df.protein_fg.values, p0=(1, 1e-6, 1))

##############################################
##############################################
# Now plot!
##############################################
##############################################

fig, ax = plt.subplots(1, 2, figsize=(5, 2.5))

# plot Schmidt and Li values
##############################################
# Load the original dataset with aboslute measurements

# Load the original dataset with aboslute measurements
data_orig = pd.read_csv('../../../data/compiled_datasets.csv')
data_orig = data_orig[data_orig.dataset == 'li_2014']

for d, df in data_orig.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
    mass = df.reported_fg_per_cell.sum()

    ax[0].plot(df.growth_rate_hr.unique(), mass, 'o', ms=4, color=colors['purple'],
            markeredgewidth=0.5, markeredgecolor='k', label=d[1])


data_orig = pd.read_csv('../../../data/compiled_datasets.csv')
data_orig = data_orig[data_orig.dataset == 'schmidt_2016']

for d, df in data_orig.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
    mass = df.reported_fg_per_cell.sum()

    ax[0].plot(df.growth_rate_hr.unique(), mass, 'o', ms=4, color=colors['light_blue'],
            markeredgewidth=0.5, markeredgecolor='k', label=d[1])

handles, labels = ax[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax[0].legend(by_label.values(), by_label.keys(), loc = 'upper left', fontsize = 6)




# approach 1 - use cell volumes from Si et al. 2017; maintain constant cellular protein concentation
##############################################
vol_glu = prot.size.lambda2size(0.58)
mass_glu = data_orig[data_orig.condition == 'glucose'].reported_fg_per_cell.sum()

for d, df in data_orig.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
    vol = prot.size.lambda2size(d[3])
    ax[1].plot(df.growth_rate_hr.unique(), (mass_glu/vol_glu) * vol , 'o', ms=4, color=colors['light_red'],
            markeredgewidth=0.5, markeredgecolor='k', label=d[1])

# plot smooth curve with data
x = np.linspace(0,2,100)
vol = prot.size.lambda2size(x)
ax[1].plot(x, (mass_glu/vol_glu) * vol, '-', alpha = 0.6,
            lw = 0.5, color = colors['light_red'], zorder = 0)

# approach 2 - use cell volumes from Si et al. 2017; maintain constant cellular protein concentation
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


schmidt_gr = data_orig.sort_values(by='growth_rate_hr', ascending = True).growth_rate_hr.unique()
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

# Plot predictions of protein along with 'data points' for Schmidt
x = np.linspace(0,2,100)
vol = prot.size.lambda2size(schmidt_gr)
ax[1].plot(schmidt_gr, pred_proteinmass_Si, 'o', ms=4, color=colors['light_green'],
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
# estimate total protein mass per cell
pred_proteinmass_Si =  (pred_RNA_protein_mass_Si/ (1+pred_RNA_protein_ratio))

#  fg per cell
ax[1].plot(x, pred_proteinmass_Si, '-', alpha = 0.6,
            lw = 0.5, color = colors['light_green'], zorder = 0)

# repeat calculations for smooth curve, predict RNA/Protein ratio >= 0.7 hr-1:
x = np.linspace(0.7, 2.0, 100)
vol = prot.size.lambda2size(x)
pred_drymass_Si = ((1.1*0.30)*(vol*1E-12)*1E15)
pred_dnamass_Si =  func(x, *popt_dna)
# Calculate total RNA + protein mass, assume DNA + RNA + protein == 90% dry mass
pred_RNA_protein_mass_Si = pred_drymass_Si*0.90 - pred_dnamass_Si
#  estimate R:P ratio
pred_RNA_protein_ratio = slope_dia_RP_B*x + intercept_dia_RP_B
# estimate total protein mass per cell
pred_proteinmass_Si =  (pred_RNA_protein_mass_Si/ (1+pred_RNA_protein_ratio))
#  fg per cell
ax[1].plot(x, pred_proteinmass_Si, '-', alpha = 0.8,
            lw = 0.9, color = colors['light_green'], zorder = 0)


# approach 3 - use measurements of protein per cell from Basan 2015
##############################################
x = np.linspace(0,2,100)
ax[1].plot(basan_df.growth_rate_hr.values, basan_df.protein_fg.values, 'o', ms=4, color=colors['light_purple'],
        markeredgewidth=0.5, markeredgecolor='k', label = 'Approach 3')
ax[1].plot(x, func(x, *popt_fg), '-', alpha = 0.6,
            lw = 0.5, color=colors['light_purple'], zorder = 0)




for ax_ in ax:
    ax_.set_ylabel('protein mass per cell [fg]', fontsize=6)
    # ax.legend(loc = 'upper left', fontsize=6)

    ax_.xaxis.set_tick_params(labelsize=5)
    ax_.yaxis.set_tick_params(labelsize=5)
    ax_.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    ax_.set_xlim([0, 2])
    ax_.set_ylim([0, 1300])


plt.tight_layout()
fig.savefig('../../../figures/schmidt_corrections_approaches.pdf', bbox_inches='tight')
