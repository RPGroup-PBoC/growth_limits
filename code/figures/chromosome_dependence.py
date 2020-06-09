#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
import prot.size
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()


data = pd.read_csv('../../data/compiled_estimate_categories.csv')
dntp = data[data['shorthand']=='dntp']
dnap = data[data['shorthand']=='dnap']
rnap = data[data['shorthand']=='rnap']
sig = data[data['shorthand']=='sigma70']
all_sig = data[data['shorthand']=='all_sigma']
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}



# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])

# Define the estimate for rRNA syntehsis
L_genes = 4500
L_spacing = 80
t_div = (2.5E6 / 600) / (60 * 60)
rate_range = np.linspace(0, 2, 200)

# Define the estimate for mRNA synthesis
t_double = np.log(2)*rate_range**-1 * 60 * 60
L_GENOME = 5E6
DENSITY = 1.1 # pg / cubic micron 
PROT_FRAC = 0.15 # Fraction of dry mass that is protein
AVG_PROT_MASS = 30E3 / (6E23 * 1E-12) # average protein mass in pg
AVG_PROT_LENGTH = 300 # in Amino Acids
AVG_AA_MASS = 110 / (6E23 * 1E-12)
volume_um3 = prot.size.lambda2size(rate_range)
N_PROT = DENSITY * volume_um3 * PROT_FRAC / (AVG_PROT_MASS)
THETA_MRNA = 1E3
R_DEG = 1 / 500
L_MRNA = 1000
R_PROD = N_MRNA * R_DEG * L_MRNA # in nt / s
R_TXN = 40 # innt / sec
R_INIT = 1
L_RNAP = 40
prediction_rRNA = 7 * L_genes * 2**((L_GENOME/2/600) / t_double) / (L_RNAP + (R_TXN / R_INIT))
prediction_mRNA =  (DENSITY * volume_um3 * PROT_FRAC * R_DEG * L_MRNA)/(AVG_PROT_MASS * R_TXN * THETA_MRNA)
prediction = prediction_rRNA + prediction_mRNA
ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in rnap.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'].values, 'o', ms=4, color=dataset_colors[g[0]],
      label=g[1], markeredgewidth=0.75, markeredgecolor='k')



ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('RNAP core enzyme per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/RNAP_trend.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])


ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in sig.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('RpoD ($\sigma_{70}$) per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/sig70_trend.pdf', bbox_inches='tight')


#%%
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
# ax.set_yscale('log')
ax.set_ylim([5E0, 15])
ax.set_xlim([0, 2])


# ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in rnap.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['concentration_uM'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('RNAP concentration [µM]', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/RNAP_conc.pdf', bbox_inches='tight')





#%%
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
ax.set_ylim([1E1, 3E4])
ax.set_xlim([0, 2])


# All sig trend

# ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in all_sig.groupby(['dataset', 'dataset_name']):
    _d = sig[(sig['dataset']==g[0]) & (sig['dataset_name']==g[1])]
    ax.plot(_d['growth_rate_hr'], _d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')
    ax.plot(d['growth_rate_hr'], d['n_complex'].values - _d['n_complex'].values, 's', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k', alpha=0.75)


ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('all $\sigma$ factors per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/all_sig_trend.pdf', bbox_inches='tight')


# %%
# predicted_dnap
prediction  = (1E7 * 2**(t_div * rate_range)) /(10 * 60 * 60 * rate_range**-1)
prediction[prediction < 1] = 1
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
# ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])


ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in dntp.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribonucleoreductases per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/rnr_trend.pdf', bbox_inches='tight')



# %%
prediction  = (5E6 * 2**(t_div * rate_range)) / (60 * 60 * 600 * rate_range**-1)
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_yscale('log')
# ax.set_ylim([1E2, 3E4])
ax.set_xlim([0, 2])


ax.plot(rate_range, prediction, 'k--', lw=0.75, label='predicted trend')
for g, d in dnap.groupby(['dataset', 'dataset_name']):
    ax.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
           label=g[1], markeredgewidth=0.75, markeredgecolor='k')

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribonucleoreductases per cell', fontsize=6)
ax.legend(fontsize=6)
plt.savefig('../../figures/rnr_trend.pdf', bbox_inches='tight')



# %%
