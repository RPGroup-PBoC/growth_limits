# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size
import prot.estimate
constants = prot.estimate.load_constants()
colors = prot.viz.plotting_style()
dataset_colors = prot.viz.dataset_colors()

data = pd.read_csv('../../data/compiled_estimate_categories.csv')

# Define constants. 
RHO = constants['density']['value']
DRY_FRAC = constants['dry_mass_frac']['value']
PROT_FRAC = DRY_FRAC * constants['theta_prot']['value']
CARB_FRAC = DRY_FRAC * constants['theta_C']['value']
VOL = constants['volume']['value']
T_DOUBLE = constants['t_double']['value']
MASS_CARB = 12/6E11 # in pg
GROWTH_RATE = constants['growth_rate']

# Define transport constants.
R_GLUC = 200 # in sugar per sec
R_XYL = 100 # in xylose per sec
R_FRUC = 200 
R_GLYC = 2000
N_GLUC = 6 # Carbons per sugar
N_XYL = 5
N_FRUC = 6
N_GLYC = 3 
# %%
# GLUCOSE_TPORTERS
N_gluc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_CARB * N_CARB * MASS_CARB * T_DOUBLE)
N_glyc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_GLYC * N_GLYC * MASS_CARB * T_DOUBLE)
N_xyl_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_XYL * N_XYL * MASS_CARB * T_DOUBLE)
N_fruc_tport = (RHO * VOL * DRY_FRAC * CARB_FRAC)\
         / (R_FRUC * N_FRUC * MASS_CARB * T_DOUBLE)


# Induced expression 
data = pd.read_csv('../../data/compiled_estimate_categories.csv', comment='#')
dataset_markers = {'li_2014':'o', 'schmidt_2016':'X',
                   'peebo_2015':'d', 'valgepea_2013':'^'}
cats = ['glucose_tport', 'glycerol_tport', 'fructose_tport', 'xylose_tport']


fig, ax = plt.subplots(2, 2, figsize=(3.5, 4))
axes = {c:a for c, a in zip(cats, ax.ravel())}
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_yscale('log')
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_ylabel('complexes per cell', fontsize=6)
    a.set_xlim([0, 2])

_ax = ax.ravel()
_ax[0].set_ylim([1E1,  1E5])
_ax[1].set_ylim([1,  5E3])
_ax[2].set_ylim([10,  5E4])
_ax[3].set_ylim([10,  5E4])

# Add the correct titles
titles = ['glucose transporters (PtsG + ManXYZ)',
          'glycerol facilitator (GlpF)',
          'fructose transporter (FruAB)',
          'xylose transporter (XylE + XylFGH)']
_ax[0].plot(GROWTH_RATE, N_gluc_tport, color='grey', lw=3, alpha=0.25,
            label='estimate')
_ax[1].plot(GROWTH_RATE, N_glyc_tport, color='grey', lw=3, alpha=0.25,
            label='estimate')
_ax[2].plot(GROWTH_RATE, N_fruc_tport, color='grey', lw=3, alpha=0.25,
            label='estimate')
_ax[3].plot(GROWTH_RATE, N_xyl_tport, color='grey', lw=3, alpha=0.25,
            label='estimate')
for a, t in zip(ax.ravel(), titles):
    prot.viz.titlebox(a, t, size=5.5, color='k', bgcolor=colors['pale_yellow'],
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
        color=_color, markeredgewidth=0.5, markeredgecolor='k', 
        label='__nolegend__')

# Add a legend. 
for g, d in data.groupby(['dataset', 'dataset_name']):
    _ax[0].plot([], [], linestyle='none', marker=dataset_markers[g[0]], color=colors['dark_green'],
                alpha=0.5, markeredgecolor='k', markeredgewidth=0.5, label=g[1],
                ms=4)

_ax[1].plot([], [], 's', color=colors['red'], markeredgecolor='k', markeredgewidth=0.5, 
            label='induced expression', ms=4)
_ax[0].legend(fontsize=5.5, loc='center')
_ax[1].legend(fontsize=5.5, loc='center')

plt.savefig('../../figures/fig2b_induced_expression.svg')

# %%
