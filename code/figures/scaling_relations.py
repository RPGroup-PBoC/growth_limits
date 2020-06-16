#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

# Define the color palettes. 
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
data = pd.read_csv('../../data/compiled_estimate_categories.csv')


# Define constants. 
RHO = 1.1 # pg / um^3
DRY_FRAC = 0.3
PROT_FRAC = 0.15
CARB_FRAC = 0.5 # Fraction of dry mass that is carbon
PHOS_FRAC = 0.03 
SULF_FRAC = 0.01
MASS_CARB = 12/6E11 # in pg
MASS_PHOS = 30/6E11 # in pg
MASS_SULF = 32/6E11 # in pg
GROWTH_RATE = np.linspace(0, 2, 200)
T_DOUBLE = 3600 * np.log(2) / GROWTH_RATE
VOL = prot.size.lambda2size(GROWTH_RATE)

# Define transport constants.
R_CARB = 200 # in sugar per sec
R_PHOS = 300 # In phosphate per sec
R_SULF = 10 # In sulfate per second
N_CARB = 6 # Carbons per sugar
N_PHOS = 1
N_SULF = 1
# %%
# GLUCOSE_TPORTERS
N_gluc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_CARB * N_CARB * MASS_CARB * T_DOUBLE)
N_phos_tport = (RHO * VOL * DRY_FRAC * PHOS_FRAC)\
         / (R_PHOS * N_PHOS * MASS_PHOS * T_DOUBLE)
N_sulf_tport = (RHO * VOL * DRY_FRAC * SULF_FRAC)\
         / (R_SULF * N_SULF * MASS_SULF * T_DOUBLE)

fig, ax = plt.subplots(3, 1, figsize=(4, 6), sharex=True)
titles = ['glucose transporters', 'phosphate transporters', 'sulfate transporters']
for i, a in enumerate(ax):
    a.set_yscale('log')
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_ylabel('complexes per cell', fontsize=6)
    prot.viz.titlebox(a, titles[i], size=6, color='k', bgcolor=colors['pale_yellow'], pad=0.02)
ax[0].set_ylim([1E2, 1E5])
ax[1].set_ylim([1E1, 5E3])
ax[2].set_ylim([1E1, 5E4])

ax[0].plot(GROWTH_RATE, N_gluc_tport, 'k--', label='prediction', lw=0.75)
ax[1].plot(GROWTH_RATE, N_phos_tport, 'k--', label='prediction', lw=0.75)
ax[2].plot(GROWTH_RATE, N_sulf_tport, 'k--', label='prediction', lw=0.75)
ax[-1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

for s, _ax in zip(['glucose_tport', 'phosphate_tport', 'sulfur_tport'], ax):
    _d = data[data['shorthand']==s]
    for g, d in _d.groupby(['dataset', 'dataset_name']):
        _ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4.5, 
                 markeredgewidth=0.5, markeredgecolor='k', 
                 color=dataset_colors[g[0]], label=g[1])

ax[0].legend(fontsize=6)
plt.savefig('../../figures/tport_system_scaling.pdf')
# %%
R_ATP = 300 # per second
C_ATP = 4
H_ATP = 4
R_PROTON = 1500
MASS_AA = 110 / 6E11 # in pg
BUDGET_FRAC = 0.8
N_synthases = (RHO * VOL * PROT_FRAC * C_ATP) / (MASS_AA * R_ATP * T_DOUBLE * BUDGET_FRAC)
N_etc = (RHO * VOL * PROT_FRAC * C_ATP * H_ATP) / (MASS_AA * R_PROTON * T_DOUBLE * BUDGET_FRAC)


fig, ax = plt.subplots(2, 1, figsize=(4, 4), sharex=True)
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_ylabel('complexes per cell', fontsize=6)
    a.set_yscale('log')

prot.viz.titlebox(ax[0], 'F1-F0 ATP synthases', size=6, bgcolor=colors['pale_yellow'], color='k')
prot.viz.titlebox(ax[1], 'electron transport chains', size=6, bgcolor=colors['pale_yellow'], color='k')

_atp = data[data['shorthand']=='atp_synthase']
_etc = data[data['shorthand']=='proton_gradient']
for g, d in _atp.groupby(['dataset', 'dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4.5, 
             markeredgewidth=0.5, markeredgecolor='k', 
             color=dataset_colors[g[0]], label=g[1])

for g, d in _etc.groupby(['dataset', 'dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4.5, 
             markeredgewidth=0.5, markeredgecolor='k', 
             color=dataset_colors[g[0]], label=g[1])

ax[-1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].plot(GROWTH_RATE, N_synthases, 'k--', lw=0.75, label='prediction')
ax[1].plot(GROWTH_RATE, N_etc, 'k--', lw=0.75, label='prediction')

plt.savefig('../../figures/energy_production_scaling.pdf')
# %%
