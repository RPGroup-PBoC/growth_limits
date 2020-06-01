#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prot.viz
colors, palette = prot.viz.bokeh_theme()
dataset_colors = {'li_2014':colors['purple'], 'schmidt_2016':colors['light_blue'],
                   'peebo_2015':colors['green'], 'valgepea_2013':colors['red']}
prot.viz.plotting_style()


def rod_SA(l,w, V):
    asp_ratio = l/w
    gamma = asp_ratio * np.pi * (asp_ratio * np.pi /4 - np.pi/12)**(-2/3)
    return gamma * V**(2/3)

def func_size(x):
    a, c = 0.53319063, -1.03724839
    return a*np.exp(-c*x)

def func_length(x):
    a, c, d = 0.49656209, -1.09303027,  1.75967254
    return a*np.exp(-c*x)+d

def func_width(x):
    a, c = 0.63830175, -0.24341639
    return a*np.exp(-c*x)

def lambda2SV(x):
    V_data = func_size(x)
    l = func_length(x)
    w = func_width(x)
    SA_data = rod_SA(l, w, V_data)

    return SA_data / V_data

def SV2lambda(x):
    a, c, d = 14.57246783,  0.38603366, -0.80393099
    return a*np.exp(-c*x)+d


# Parameters and calculations #
###############################

V = np.linspace(0.5,50, 500)

SA_rod = 2 * np.pi *  V**(2/3)
SA_V_ratio_rod = SA_rod / V

# Sphere V = (4/3) pi r**3
# Sphere SA = 4 pi r**2
SA_sphere = V**(2/3) * ((4/3) * np.pi)**(-2/3) * 4 * np.pi
SA_V_ratio_sphere = SA_sphere / V

# ATP equivalents demand w.r.t. volume ; 1E6 ATP/(um3 s)
Pv = 1E6 * V

# ATP max - half surface area devoted to respiration
Ps_um_resp = ((3)/ (1E-6))
Ps_resp_rod = Ps_um_resp * SA_rod * 0.5
Ps_resp_sphere = Ps_um_resp * SA_sphere * 0.5

# # ATP max - half surface area devoted to fermentation
# Ps_um_ferm= ((180*2)/ (50E-6))
# Ps_ferm_rod = Ps_um_ferm * SA_rod * 0.5
# Ps_ferm_sphere = Ps_um_ferm * SA_sphere * 0.5

#### for Fill_between, need common x-axis array
SA_ = np.linspace(np.min(SA_sphere), np.max(SA_rod), 500)
V_sphere = (SA_/(((4/3) * np.pi)**(-2/3) * 4 * np.pi))**(3/2)
SA_V_ratio_sphere_ = SA_/V_sphere

V_rod = (SA_/(2 * np.pi))**(3/2)
SA_V_ratio_rod_ = SA_/V_rod

# ATP max - half surface area devoted to respiration
Ps_um_resp = ((3)/ (1E-6))
Ps_resp_ = Ps_um_resp * SA_ * 0.5

###############################
###############################

# %%
fig, ax = plt.subplots(1, 2, figsize=(3, 2),
        gridspec_kw={'width_ratios': [6, 1]})#}, sharey='all')

# plot the demand
ax[0].plot(Pv, SA_V_ratio_rod, color=colors['dark_green'], label='rod',
           alpha=0.9, lw = 0.5, ls = '-.')
ax[0].plot(Pv, SA_V_ratio_sphere, color=colors['dark_green'], label='sphere',
           alpha=0.9, lw = 0.5, ls = '--')

ax[0].fill_between(Pv, y1 = SA_V_ratio_sphere, y2=SA_V_ratio_rod,
        color=colors['dark_green'],alpha=0.2, lw = 0)

# plot the max for respiration
ax[0].plot(Ps_resp_rod, SA_V_ratio_rod, color=colors['blue'],
            label='rod', alpha=0.9, lw = 0.5, ls = '-.')
ax[0].plot(Ps_resp_sphere, SA_V_ratio_sphere, color=colors['blue'],
            label='sphere', alpha=0.9, lw = 0.5, ls = '--')

ax[0].fill_between(Ps_resp_, y1 = SA_V_ratio_sphere_, y2 = SA_V_ratio_rod_,
        color=colors['blue'],alpha=0.2, lw = 0)

# # plot the max for fermentation
# ax[0].plot(Ps_ferm_rod, SA_V_ratio_rod, color=colors['red'],
#             label='rod', alpha=0.9)
# ax[0].plot(Ps_ferm_sphere, SA_V_ratio_sphere, color=colors['red'],
#             label='sphere', alpha=0.4)

#%
# # Populate second plot with growth rates
# S/V for E. coli datasets
# Load the data set
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

for g, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    SV = lambda2SV(g[2])
    ax[1].plot(1, SV, 'o', color=dataset_colors[g[0]],
                    alpha=0.75, markeredgecolor='k', markeredgewidth=0.25,
                    label = g[2], ms=4, zorder=10)

# Format the axes
for a in ax:
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
    a.set_xscale('log')
    a.set_ylim([1.5, 8.0])
    # a.legend(fontsize=5, loc='lower right')

ax[0].set_xlim([np.min(Pv), np.max(Ps_resp_)])
ax[0].set_xlabel('ATP equivalents per s', fontsize=6)
ax[0].set_ylabel('S/V ratio [$\mu$m]', fontsize=6)

ax[1].xaxis.set_ticks([])
ax[1].set_yticks(lambda2SV(np.array([0,0.5,1,2])))
ax[1].set_yticklabels(np.array([0,0.5,1,2]))
ax[1].yaxis.set_label_position("right")
ax[1].yaxis.tick_right()
ax[1].set_ylabel('measured growth rate [hr$^{-1}$]', fontsize=6)
ax[1].xaxis.set_tick_params(labelsize=5)
ax[1].yaxis.set_tick_params(labelsize=5)

# move second plot closer
box = ax[1].get_position()
box.x0 = box.x0 + 0.045
box.x1 = box.x1 + 0.045
ax[1].set_position(box)

plt.tight_layout()
plt.savefig('../../figures/energy_estimate__SV_scaling_plots.pdf')
# %%
