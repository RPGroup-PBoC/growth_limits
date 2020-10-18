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
R_XYL = 100 # in xylose per sec
R_FRUC = 200 
R_GLYC = 2000
R_PHOS = 300 # In phosphate per sec
R_SULF = 10 # In sulfate per second
N_CARB = 6 # Carbons per sugar
N_XYL = 5
N_FRUC = 6
N_GLYC = 3 
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
N_glyc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_GLYC * N_GLYC * MASS_CARB * T_DOUBLE)
N_xyl_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_XYL * N_XYL * MASS_CARB * T_DOUBLE)
N_fruc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_FRUC * N_FRUC * MASS_CARB * T_DOUBLE)



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
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_yscale('log')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('number of phosphate transporters', fontsize=6)
ax.set_xlim([0, 2])

ax.plot(GROWTH_RATE, N_phos_tport, 'gray', lw=5, alpha=0.25, label='cell size dependence')
ax.plot(0.5, 133, 'o', color=colors['dark_brown'], alpha=0.4, ms=6, label='point estimate', zorder=1000)
ax.hlines(133, 0, 0.5, 'k', linestyle='--', lw=0.5, label='__nolegend__', zorder=999)
ax.vlines(0.5, 0, 133, 'k', linestyle='--', lw=0.5, label='__nolegend__', zorder=999)

for g, d in data[data['shorthand']=='phosphate_tport'].groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', markeredgewidth=0.5,
    markeredgecolor='k', markerfacecolor=dataset_colors[g[0]], label=g[1],
    ms=4, alpha=0.75)

ax.legend(fontsize=6)
plt.savefig('../../figures/phospho_tport_scaling_plot.pdf')

# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_yscale('log')
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('number of sulfate transporters', fontsize=6)
ax.set_xlim([0, 2])
ax.set_ylim([1E1, 1E5])

ax.plot(GROWTH_RATE, N_sulf_tport, 'gray', lw=5, alpha=0.25, label='cell size dependence')
ax.plot(0.5, 1E3, 'o', color=colors['dark_brown'], alpha=0.4, ms=6, label='point estimate', zorder=1000)
ax.hlines(1E3, 0, 0.5, 'k', linestyle='--', lw=0.5, label='__nolegend__', zorder=999)
ax.vlines(0.5, 0, 1E3, 'k', linestyle='--', lw=0.5, label='__nolegend__', zorder=999)

for g, d in data[data['shorthand']=='sulfur_tport'].groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', markeredgewidth=0.5,
    markeredgecolor='k', markerfacecolor=dataset_colors[g[0]], label=g[1],
    ms=4, alpha=0.75)

ax.legend(fontsize=6, ncol=2, handlelength=1)
plt.savefig('../../sulfur_tport_scaling_plot.svg')

# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.set_ylim([1E1, 1E4])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('phosphate transporters per cell', fontsize=6)
ax.plot(0.5, 133, 'o', color=colors['dark_brown'], alpha=0.5, ms=6, label='point estimate', zorder=100)
ax.hlines(133, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__', zorder=99)
ax.vlines(0.5, 0, 133, 'k', linestyle='--', lw=0.75, label='__nolegend__', zorder=99)
for g, d in data[data['shorthand']=='phosphate_tport'].groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4.5, markeredgecolor='k', 
            markeredgewidth=0.5, color=dataset_colors[g[0]], label=g[1])

ax.legend(fontsize=6)
plt.savefig('../../figures/phosphate_point_estimate.pdf', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_xlim([0, 2])
ax.set_ylim([1E1, 1E4])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('phosphate transporters per cell', fontsize=6)
ax.plot(GROWTH_RATE, N_phos_tport, 'k--',  label='estimate', zorder=100)
for g, d in data[data['shorthand']=='phosphate_tport'].groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4.5, markeredgecolor='k', 
            markeredgewidth=0.5, color=dataset_colors[g[0]], label=g[1])
ax.legend(fontsize=6)
plt.savefig('../../figures/phosphate_continuum_estimate.pdf', bbox_inches='tight')


# %%
# Induced expression 
data = pd.read_csv('../../data/compiled_estimate_categories.csv', comment='#')
dataset_markers = {'li_2014':'o', 'schmidt_2016':'X',
                   'peebo_2015':'d', 'valgepea_2013':'^'}
cats = ['carbon_tport_tot', 'glucose_tport', 'glycerol_tport', 'fructose_tport', 'xylose_tport']


fig, ax = plt.subplots(2, 3, figsize=(6.5, 4))
axes = {c:a for c, a in zip(cats, ax.ravel())}
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_yscale('log')
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_ylabel('complexes per cell', fontsize=6)
    a.set_xlim([0, 2])

_ax = ax.ravel()
_ax[-1].axis('off')
_ax[0].set_ylim([5E2,  1E6])
_ax[1].set_ylim([5E2,  1E6])
_ax[2].set_ylim([1,  5E3])
_ax[3].set_ylim([10,  1E4])
_ax[4].set_ylim([10,  1E4])

# Add the correct titles
titles = ['all carbohydrate transporters',
          'glucose transporters (PtsG + ManXYZ)',
          'glycerol facilitator (GlpF)',
          'fructose transporter (FruAB)',
          'xylose transporter (XylE + XylFGH)']
_ax[1].plot(GROWTH_RATE, N_gluc_tport, color='grey', lw=3, alpha=0.25)
_ax[2].plot(GROWTH_RATE, N_glyc_tport, color='grey', lw=3, alpha=0.25)
_ax[3].plot(GROWTH_RATE, N_fruc_tport, color='grey', lw=3, alpha=0.25)
_ax[4].plot(GROWTH_RATE, N_xyl_tport, color='grey', lw=3, alpha=0.25)
for a, t in zip(ax.ravel(), titles):
    prot.viz.titlebox(a, t, size=6, color='k', bgcolor=colors['pale_yellow'],
                    boxsize=0.12)



for g, d in data.groupby(['dataset', 'growth_rate_hr', 'condition']):
    for c, a in axes.items():
        _d = d[d['shorthand']==c]
        if (c.split('_tport')[0] in g[-1]) & ('glucose' not in g[-1]):
            _color = colors['red']
            alpha=1

            if 'pAA' in g[-1]:
                label = 'glycerol + A.A.'
            else:
                label = c.split('_tport')[0]
            a.text(_d['growth_rate_hr'] + 0.05, _d['n_complex'],
            label, fontsize=6)
        else:
            alpha = 0.5
            _color = colors['dark_green'] 

        marker = dataset_markers[g[0]]
        a.plot(_d['growth_rate_hr'], _d['n_complex'], linestyle='none',
              marker=marker, ms=4, alpha=alpha,
        color=_color, markeredgewidth=0.5, markeredgecolor='k')

# Add a legend. 
for g, d in data.groupby(['dataset', 'dataset_name']):
    _ax[-1].plot([], [], linestyle='none', marker=dataset_markers[g[0]], color=colors['dark_green'],
                alpha=0.5, markeredgecolor='k', markeredgewidth=0.5, label=g[1],
                ms=4)
_ax[-1].plot([], [], 's', color=colors['red'], markeredgecolor='k', markeredgewidth=0.5, 
            label='induced expression', ms=4)
_ax[-1].legend(fontsize=7.5, loc='center')
plt.tight_layout()
plt.savefig('../../figures/scaling_induced_expression.pdf')
# %%
