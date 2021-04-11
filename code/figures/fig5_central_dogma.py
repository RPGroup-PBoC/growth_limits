
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

# Set up the figure canvas.
# DNA replication,  dNTP,
# RNA synthesis  legend
# protein synthesis (2),
# additional column:
# DNA pol conc. , rna pol RpoD, ribo  synthesis time

# ribosomal synthesis time schematic

fig = plt.figure(figsize = (8,6))
widths = [2, 2, 2]
heights = [2, 2, 2]
spec = fig.add_gridspec(ncols=3, nrows=3, width_ratios=widths,
                          height_ratios=heights)

ax1 = fig.add_subplot(spec[0, 0]) # DNA pol III
ax2 = fig.add_subplot(spec[0, 1]) # dNTP
ax3 = fig.add_subplot(spec[1, 0]) # RNA pol
ax4 = fig.add_subplot(spec[2, 0]) # ribosomes
ax5 = fig.add_subplot(spec[2, 1]) # tRNA

ax6 = fig.add_subplot(spec[0, 2]) # DNA pol conc.
ax7 = fig.add_subplot(spec[1, 2]) # RpoD

###################################
# dNTP
###################################
dnap = data[data['shorthand']=='dnap']
rnr = data[data['shorthand']=='dntp']

# Compute the cell size dependence.
growth_rate = constants['growth_rate']['value']
t_double = constants['t_double']['value']
n_ori = constants['N_ori']['value']
L_genome = 4.6E6 # in nt
r_rnr = 10
r_dna = 600 # in nt per sec
n_pol = 2 # per replication fork
n_fork = 2
N_rnr = 2 * n_ori * L_genome  / (r_rnr * t_double)
N_dnap = n_fork * n_pol * n_ori

ax2.xaxis.set_tick_params(labelsize=8)
ax2.yaxis.set_tick_params(labelsize=8)
ax2.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
ax2.set_ylabel('ribonucleotide reductases\nper cell', fontsize=8)

# Format the axes
ax2.set_yscale('log')
ax2.set_ylim([1E1, 1E4])
ax2.set_xlim([0, 2])

ax2.plot(growth_rate[growth_rate > 0.23], N_rnr[growth_rate > 0.23],'-', color='grey', lw=3, alpha=0.5, label='replication fork dependence')
ax2.plot(growth_rate[growth_rate <= 0.23], N_rnr[growth_rate <= 0.23], ':', color='grey', lw=3, alpha=0.5, label='__nolegend__')


# Plot the predictions
ax2.plot(0.5, 200, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax2.hlines(250, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax2.vlines(0.5, 10, 200, 'k', linestyle='--', lw=0.75, label='__nolegend__')

for g, d in rnr.groupby(['dataset', 'dataset_name']):
    ax2.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])


###################################
# DNA replication
###################################

for a in [ax1,ax6]:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
ax1.set_ylabel('DNA polymerase III\nper cell', fontsize=8)
ax6.set_ylabel('DNA polymerase III\nconcentration [nM]', fontsize=8)

# Format the axes
ax1.set_yscale('log')
ax6.set_yscale('log')
ax1.set_ylim([1, 5E2])
ax1.set_xlim([0, 2])
ax6.set_ylim([1, 500])
ax6.set_xlim([0, 2])


# Plot the predictions
ax1.plot(growth_rate, N_dnap, '-', lw=3, color='grey', alpha=0.3, label='replication fork dependence')
ax1.plot(growth_rate[growth_rate > 0.23], N_dnap[growth_rate > 0.23], '-', lw=3, color='grey', alpha=0.3, label='replication fork dependence')
ax1.plot(growth_rate[growth_rate <= 0.23], N_dnap[growth_rate <= 0.23], ':', lw=3, color='grey', alpha=0.3, label='__nolegend__')
ax1.plot(0.5, 6.5, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='estimated value')
ax1.hlines(6, 0, 0.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')
ax1.vlines(0.5, 1, 6.5, 'k', linestyle='--', lw=0.75, label='__nolegend__')



# Plot the concentration range of 50 - 200 nM
ax6.fill_between([0, 2], 50, 200, color='k', alpha=0.25, label='DNA Pol. III H.E. binding affinity\n (Ason et al. 2000)')

for g, d in dnap.groupby(['dataset', 'dataset_name']):
    ax1.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label=g[1])
    ax6.plot(d['growth_rate_hr'], d['concentration_uM'] * 1E3, 'o', ms=4, color=dataset_colors[g[0]],
            alpha=0.75, markeredgewidth=0.5, markeredgecolor='k', label='__nolegend__')
# ax6.legend(fontsize=8)
# ax1.legend(fontsize=8, ncol=2)


###################################
# RNA synthesis
###################################

rnap = data[data['shorthand']=='rnap']
sig70 = data[data['shorthand']=='sigma70']


# Compute the predicted trends..
growth_rate = constants['growth_rate']['value']
N_ori = constants['N_ori']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
t_double = constants['t_double']['value']

m_aa = 110 / 6E11 # in pg
m_prot = (300 * 110) / 6E11 # in pg
n_prot = (theta_prot * theta_dry * cell_mass) /m_prot
gamma_mRNA = 1 / 300
n_mRNA = (n_prot / 1E3) * gamma_mRNA
L_mRNA = 1000
L_rRNA = 4500
r_txn = 40 # in nt/s
n_operon = 7
footprint = 80 # in nt
n_tRNA = theta_prot * theta_dry * cell_mass / (m_aa * t_double)
L_tRNA = 80

N_polymerase = (L_rRNA * n_operon  * N_ori/ footprint) + (n_mRNA * L_mRNA / r_txn) + (n_tRNA * L_tRNA) / (r_txn * t_double)


for a in [ax3, ax7]:
    a.plot(growth_rate[growth_rate > 0.23], N_polymerase[growth_rate > 0.23], '-', lw=3, color='grey', alpha=0.4, label='replication fork scaling')
    a.plot(growth_rate[growth_rate <= 0.23], N_polymerase[growth_rate <= 0.23], ':', lw=3, color='grey', alpha=0.4, label='__nolegend__')
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
    a.set_yscale('log')
    a.set_xlim([0, 2])
ax3.set_ylabel('RNA Polymerases\nper cell', fontsize=8)
ax7.set_ylabel('$\sigma^{70}$ (RpoD)\nper cell', fontsize=8)
ax3.set_yticks([1E2, 1E3, 1E4, 1E5])
ax3.set_ylim([1E2, 1E5])
ax7.set_yticks([1E2, 1E3, 1E4])
ax7.set_ylim([1E2, 1E4])

# Plot the predictions
for a in [ax3, ax7]:
    a.plot(0.5, 1E3, 'o', ms=6, alpha=0.4, color=colors['dark_brown'], label='estimated value')
    a.vlines(0.5, 1, 1E3, color='k', linestyle='--', label='__nolegend__', lw=0.75)
    a.hlines(1E3, 0, 0.5, color='k', linestyle='--', label='__nolegend__', lw=0.75)

# plot the data
for p, a in zip([rnap, sig70], [ax3, ax7]):
    for g, d in p.groupby(['dataset', 'dataset_name']):
        a.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
                markeredgewidth=0.5, markeredgecolor='k', label=g[1])

# Add legends.
for a in [ax3, ax7]:
    # a.legend(fontsize=8, loc='lower right')
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

###################################
# protein synthesis
###################################
data = pd.read_csv('../../data/compiled_estimate_categories.csv', comment='#')
data = data[data['shorthand']=='ribosome']

# Compute the trend
growth_rate = constants['growth_rate']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
t_double = constants['t_double']['value']
m_aa = 110 / 6E11 # in pg
k_tsl = 15 # per sec

N_ribosomes = cell_mass * theta_dry * theta_prot / (m_aa * k_tsl * t_double)

ax4.xaxis.set_tick_params(labelsize=8)
ax4.yaxis.set_tick_params(labelsize=8)
ax4.set_yscale('log')
ax4.set_xlim([0, 2])
ax4.set_ylim([1E3, 3E5])
ax4.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
ax4.set_ylabel('ribosomes\nper cell', fontsize=8)

# Plot the cell size dependence
ax4.plot(growth_rate[growth_rate > 0.23], N_ribosomes[growth_rate > 0.23], lw=3, color='grey', alpha=0.4, label='cell size dependence')
ax4.plot(growth_rate[growth_rate <= 0.23], N_ribosomes[growth_rate <= 0.23], ':', lw=3, color='grey', alpha=0.4, label='__nolegend__')

# Plot the point estimate.
estimate = 1E4
ax4.plot(0.5, estimate, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax4.vlines(0.5, 0, estimate, 'k', linestyle='--', lw=1, label='__nolegend__')
ax4.hlines(estimate, 0, 0.5, 'k', linestyle='--', lw=1, label='__nolegend__')

# Plot the experimetnal data
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax4.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])


###################################
# tRNA synthesis
###################################
data = pd.read_csv('../../data/compiled_estimate_categories.csv', comment='#')
data = data[data['shorthand']=='trna']

# Compute the trend
growth_rate = constants['growth_rate']['value']
cell_mass = constants['cell_mass']['value']
theta_dry = constants['dry_mass_frac']['value']
theta_prot = constants['theta_prot']['value']
t_double = constants['t_double']['value']
m_aa = 110 / 6E11 # in pg
k_trna = 20 # per sec

N_synthase = cell_mass * theta_dry * theta_prot / (m_aa * k_trna * t_double)

ax5.xaxis.set_tick_params(labelsize=8)
ax5.yaxis.set_tick_params(labelsize=8)
ax5.set_yscale('log')
ax5.set_xlim([0, 2])
ax5.set_ylim([1E3, 3E5])
ax5.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
ax5.set_ylabel('tRNA synthetases\nper cell', fontsize=8)

# Plot the cell size dependence
ax5.plot(growth_rate[growth_rate > 0.23], N_synthase[growth_rate > 0.23], lw=3, color='grey', alpha=0.4, label='cell size dependence')
ax5.plot(growth_rate[growth_rate <= 0.23], N_synthase[growth_rate <= 0.23], ':', lw=3, color='grey', alpha=0.4, label='__nolegend__')

# Plot the point estimate.
estimate = 1E4
ax5.plot(0.5, estimate, 'o', ms=4.5, color=colors['dark_brown'], alpha=0.4, label='point estimate')
ax5.vlines(0.5, 0, estimate, 'k', linestyle='--', lw=1, label='__nolegend__')
ax5.hlines(estimate, 0, 0.5, 'k', linestyle='--', lw=1, label='__nolegend__')

# Plot the experimetnal data
for g, d in data.groupby(['dataset', 'dataset_name']):
    ax5.plot(d['growth_rate_hr'], d['n_complex'], 'o', ms=4, color=dataset_colors[g[0]],
            markeredgewidth=0.5, markeredgecolor='k', label=g[1])



# ax2.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig('../../figures/fig4_central_dogma.pdf')#, bbox_inches='tight')





# %%
