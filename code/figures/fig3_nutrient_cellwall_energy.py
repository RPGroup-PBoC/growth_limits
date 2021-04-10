
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()
constants = prot.estimate.load_constants()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')
carbon = data[data['shorthand']=='carbon_tport']
phos = data[data['shorthand']=='phosphate_tport']
sulf = data[data['shorthand']=='sulfur_tport']
lipid = data[data['shorthand']=='fas']
pg = data[data['shorthand']=='transpeptidases']
proton = data[data['shorthand']=='proton_gradient']
atp = data[data['shorthand']=='atp_synthase']

###################################
# CARBON TRANSPORT (Fig 3A)
###################################

# Set up the figure canvas.
fig = plt.figure(figsize = (12,6))
widths = [2, 2, 2, 2]
heights = [2, 2]
spec = fig.add_gridspec(ncols=4, nrows=2, width_ratios=widths,
                          height_ratios=heights)
ax = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[0, 2])

ax8 = fig.add_subplot(spec[1, 0])
ax9 = fig.add_subplot(spec[1, 1])

ax10 = fig.add_subplot(spec[1, 2])
ax11 = fig.add_subplot(spec[1, 3])


# Define constants
theta_C = constants['dry_mass_frac']['value'] * constants['theta_C']['value']
rho = constants['density']['value']
vol = constants['volume']['value']
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
mass = constants['cell_mass']['value']
m_carbon = 12/6E11 # in pg
r_carbon = 1000 # in C / s
N_tporters = (theta_C * mass)/ (m_carbon * r_carbon * t_double)

ax.xaxis.set_tick_params(labelsize=10)
ax.yaxis.set_tick_params(labelsize=10)
ax.set_xlim([0, 2])
ax.set_ylim([1E2, 5E5])
ax.set_yscale('log')
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax.set_ylabel('carbon transporters\nper cell (PTS)', fontsize=10)

# Plot the scaling argument
ax.plot(0.5, 2E3, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax.vlines(0.5, 1E2, 2E3, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax.hlines(2E3, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in carbon.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax.plot(growth_rate[growth_rate > 0.23], N_tporters[growth_rate > 0.23], '-', lw=3, color='grey', label='cell size dependence',
alpha=0.3)
ax.plot(growth_rate[growth_rate <= 0.23], N_tporters[growth_rate <= 0.23], ':', lw=3, color='grey', label='__nolegend__',
alpha=0.3)


# data = pd.read_csv('../../data/compiled_estimate_categories.csv')
# data = data[data['shorthand']=='phosphate_tport']

# Define constants
theta_P = constants['dry_mass_frac']['value'] * constants['theta_P']['value']
rho = constants['density']['value']
vol = constants['volume']['value']
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
mass = constants['cell_mass']['value']
m_phos = 30/5E11 # in pg
r_phos = 300 # in P / s
N_tporters = (theta_P * mass)/ (m_phos * r_phos * t_double)

# data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# # Define constants.
# RHO = constants['density']['value']
# DRY_FRAC = constants['dry_mass_frac']['value']
# PROT_FRAC = DRY_FRAC * constants['theta_prot']['value']
# CARB_FRAC = DRY_FRAC * constants['theta_C']['value']
# VOL = constants['volume']['value']
# T_DOUBLE = constants['t_double']['value']
# MASS_CARB = 12/6E11 # in pg
# GROWTH_RATE = constants['growth_rate']['value']

# Define transport constants.
# R_GLUC = 200 # in sugar per sec
# R_XYL = 50 # in xylose per sec
# R_FRUC = 200
# R_GLYC = 2000
# N_GLUC = 6 # Carbons per sugar
# N_XYL = 5
# N_FRUC = 6
# N_GLYC = 3
# # GLUCOSE_TPORTERS
# N_gluc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
#          / (R_GLUC * N_GLUC * MASS_CARB * T_DOUBLE)
# N_glyc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
#          / (R_GLYC * N_GLYC * MASS_CARB * T_DOUBLE)
# N_xyl_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
#          / (R_XYL * N_XYL * MASS_CARB * T_DOUBLE)
# N_fruc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
#          / (R_FRUC * N_FRUC * MASS_CARB * T_DOUBLE)

##################################
# PHOSPHATE TRANSPORT (Fig 3B)
###################################
# data = pd.read_csv('../../data/compiled_estimate_categories.csv')
ax2.xaxis.set_tick_params(labelsize=10)
ax2.yaxis.set_tick_params(labelsize=10)
ax2.set_xlim([0, 2])
ax2.set_ylim([1E1, 1E4])
ax2.set_yscale('log')
ax2.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax2.set_ylabel('phosphate transporters\nper cell (PitA + PitB)', fontsize=10)

# Plot the scaling argument
ax2.plot(0.5, 2E2, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='point estimate',
        zorder=1000)
ax2.vlines(0.5, 1E1, 2E2, color='k', linestyle='--', lw=0.75, label='__nolegend__',
        zorder=999)
ax2.hlines(2E2, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__',
        zorder=999)

# Plot the data
for g, d in phos.groupby(['dataset', 'dataset_name']):
    ax2.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax2.plot(growth_rate[growth_rate > 0.23], N_tporters[growth_rate > 0.23], '-', lw=3, color='grey', label='cell size dependence',
alpha=0.3)
ax2.plot(growth_rate[growth_rate <= 0.23], N_tporters[growth_rate <= 0.23], ':', lw=3, color='grey', label='__nolegend__',
alpha=0.3)

# data = pd.read_csv('../../data/compiled_estimate_categories.csv')
###################################
# SULFUR TRANPSORT (Fig 3 C)
###################################


# Define constants
theta_S = constants['dry_mass_frac']['value'] * constants['theta_S']['value']
rho = constants['density']['value']
vol = constants['volume']['value']
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
mass = constants['cell_mass']['value']
m_sulf = 32/6E11 # in pg
r_sulf = 10 # in C / s
N_tporters = (theta_S * mass)/ (m_sulf * r_sulf * t_double)

ax3.xaxis.set_tick_params(labelsize=10)
ax3.yaxis.set_tick_params(labelsize=10)
ax3.set_xlim([0, 2])
ax3.set_ylim([1E1, 5E4])
ax3.set_yscale('log')
ax3.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax3.set_ylabel('sulfate transporters\nper cell (CysUWA)', fontsize=10)

# Plot the scaling argument
ax3.plot(0.5, 1E3, 'o', ms=6, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax3.vlines(0.5, 1E1, 1E3, color='k', linestyle='--', lw=0.75, label='__nolegend__')
ax3.hlines(1E3, 0, 0.5, color='k', linestyle='--', lw=0.75, label='__nolegend__')

# Plot the data
for g, d in sulf.groupby(['dataset', 'dataset_name']):
    ax3.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])

ax3.plot(growth_rate[growth_rate > 0.23], N_tporters[growth_rate > 0.23], '-', lw=3, color='grey', label='cell size dependence',
alpha=0.3)
ax3.plot(growth_rate[growth_rate <= 0.23], N_tporters[growth_rate <= 0.23], ':', lw=3, color='grey', label='__nolegend__',
alpha=0.3)

##########################################
# CELL ENVELOPEE BIOSYNTHESIS (Fig. 3 D-E)
##########################################
# data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# Compute the scaling relations.
growth_rate = constants['growth_rate']['value']
surface_area = constants['surface_area']['value']
t_double = constants['t_double']['value']

rho_pg = 0.25 # in pg per fL
A_lipid = 0.5 / 1E6 # in square microns
N_leaflet = 4
kcat_acp = 1
kcat_tpd = 2
xlink_frac = 0.2
w_pg = 0.005 # thickness of pg in um
m_murein = 1000 / 6E11 # Mass of murein monomer in pg
theta_lipid = 0.4
N_fabs = (theta_lipid * surface_area * N_leaflet) / (A_lipid * kcat_acp * t_double)
N_tpds = (xlink_frac * w_pg * surface_area * rho_pg) / (m_murein * kcat_tpd * t_double)

# Generate the figures
ax8.set_xlim([0, 2])
ax8.set_yscale('log')
ax8.set_ylim([1E2, 1E5])
ax8.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax8.set_ylabel('ACP dehydratases per cell\n(FabZ + FabA)', fontsize=10)
ax8.xaxis.set_tick_params(labelsize=10)
ax8.yaxis.set_tick_params(labelsize=10)

# Plot the scaling relationship
ax8.plot(growth_rate[growth_rate > 0.23], N_fabs[growth_rate > 0.23], '-', lw=3, color='grey', label='surface area scaling',
        alpha=0.4)
ax8.plot(growth_rate[growth_rate <= 0.23], N_fabs[growth_rate <= 0.23], ':', lw=3, color='grey', label='surface area scaling',
        alpha=0.4)

# Plot the prediction
ax8.plot(0.5, 4E3, 'o', ms=5,  color=colors['dark_brown'],  label='point estimate',
        alpha=0.4)
ax8.hlines(4E3, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax8.vlines(0.5, 0, 4E3, 'k', linestyle='--', lw=0.75, label='__nolegend__')

# plot the data
for g, d in lipid.groupby(['dataset', 'dataset_name']):
    ax8.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4,
            color=dataset_colors[g[0]], alpha=0.75, markeredgecolor='k',
            markeredgewidth=0.5, label=g[1])
###################################
# lipid
###################################

ax9.set_xlim([0, 2])
ax9.set_yscale('log')
ax9.set_ylim([5E0, 5E3])
ax9.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax9.set_ylabel('murein transpeptidases per cell', fontsize=10)
ax9.xaxis.set_tick_params(labelsize=10)
ax9.yaxis.set_tick_params(labelsize=10)

# Plot the scaling relationship
ax9.plot(growth_rate[growth_rate > 0.23], N_tpds[growth_rate > 0.23], '-', lw=3, color='grey', label='surface area scaling',
        alpha=0.4)
ax9.plot(growth_rate[growth_rate <= 0.23], N_tpds[growth_rate <= 0.23], ':', lw=3, color='grey', label='__nolegend__',
        alpha=0.4)

# Plot the prediction
ax9.plot(0.5, 100, 'o', ms=5,  color=colors['dark_brown'],  label='point estimate',
        alpha=0.4)
ax9.hlines(100, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax9.vlines(0.5, 0, 100, 'k', linestyle='--', lw=0.75, label='__nolegend__')

# plot the data
for g, d in pg.groupby(['dataset', 'dataset_name']):
    ax9.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4,
            color=dataset_colors[g[0]], alpha=0.75, markeredgecolor='k',
            markeredgewidth=0.5, label=g[1])


###################################
# ATP SYNTHESIS (Fig. 3 F)
###################################
# Compute the scaling trend.
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']

m_aa = 110 / 6E11 # in pg
r_atp = 300 # per second per synthase
atp_aa = 5

tot_prot = prot.size.lambda2P(growth_rate) / 1E3
N_synthase = (tot_prot * atp_aa) / (m_aa * r_atp * t_double / np.log(2))

ax10.xaxis.set_tick_params(labelsize=10)
ax10.yaxis.set_tick_params(labelsize=10)
ax10.set_yscale('log')
ax10.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax10.set_ylabel('ATP synthases per cell', fontsize=10)
ax10.set_xlim([0, 2])
ax10.set_ylim([1E2, 5E4])

# Plot the scaling relationship
ax10.plot(growth_rate[growth_rate > 0.23], N_synthase[growth_rate >= 0.23], '-', lw=3, color='grey', alpha=0.4, label='cell size dependence')
ax10.plot(growth_rate[growth_rate <= 0.23], N_synthase[growth_rate <= 0.23], ':', lw=3, color='grey', alpha=0.4, label='__nolegend_-')

# Plot the estimate value
ax10.plot(0.5, 3000, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.5, label='estimated value')
ax10.vlines(0.5, 0, 3000, 'k', lw=1, linestyle='--', label='__nolegend__')
ax10.hlines(3000, 0, 0.5, 'k', lw=1, linestyle='--', label='__nolegend__')

for g, d in atp.groupby(['dataset', 'dataset_name', 'condition']):
    ax10.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])

###################################
# PROTON GRADIENT (Fig. 3 G)
###################################
# data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# Compute the scaling trend.
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
m_aa = 110 / 6E11 # in pg
r_atp = 300 # per second per synthase
atp_aa = 5
prot_atp = 4
r_etc = 1500

tot_prot = prot.size.lambda2P(growth_rate) / 1E3
N_synthase = (tot_prot * atp_aa) / (m_aa * r_atp * t_double / np.log(2))
N_ETC = N_synthase * prot_atp * r_atp / r_etc

ax11.xaxis.set_tick_params(labelsize=10)
ax11.yaxis.set_tick_params(labelsize=10)
ax11.set_yscale('log')
ax11.set_xlabel('growth rate [hr$^{-1}$]', fontsize=10)
ax11.set_ylabel('electron transport complexes\nper cell', fontsize=10)
ax11.set_xlim([0, 2])
ax11.set_ylim([1E2, 5E4])
# Plot the scaling relationship
ax11.plot(growth_rate[growth_rate > 0.23], N_ETC[growth_rate > 0.23], '-', lw=3, color='grey', alpha=0.4, label='cell size dependence')
ax11.plot(growth_rate[growth_rate <= 0.23], N_ETC[growth_rate <= 0.23], ':', lw=3, color='grey', alpha=0.4, label='__nolegend__')

# Plot the estimate value
ax11.plot(0.5, 2500, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.5, label='estimated value')
ax11.vlines(0.5, 0, 2500, 'k', lw=1, linestyle='--', label='__nolegend__')
ax11.hlines(2500, 0, 0.5, 'k', lw=1, linestyle='--', label='__nolegend__')

for g, d in proton.groupby(['dataset', 'dataset_name']):
    ax11.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])

plt.tight_layout()
plt.savefig('../../figures/fig3_nutrient_cellwall_energy.pdf')