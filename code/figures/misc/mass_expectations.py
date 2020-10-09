import numpy as np
import pandas as pd
from scipy import stats
import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import prot.viz
colors, palette = prot.viz.bokeh_theme()
prot.viz.plotting_style()

density = 1.1
frac_dry = 0.3
growth_rate = np.linspace(0,2, 10)

RNA_protein_ratio_data = [0.2, 0.4]
RNA_protein_ratio_data_x = [0.5, 1.5]

slope, intercept, _, _, _ = stats.linregress(RNA_protein_ratio_data_x,
                                                RNA_protein_ratio_data)

RNA_protein_ratio = slope * growth_rate + intercept

# Biology by the numbers,  Sinauer, 1990 (BNID 104954)
# growth rate ln2/(40/60) ~ 1 ; their RNA/protein is ~0.36 measured which seems on target
# frac of mass that is protein is 55 %
# frac mass that is protein and RNA is 75 %

vol = 0.27 * 2**(1.1  * growth_rate / np.log(2))

mass_dry = vol * frac_dry * density * 1000
mass_RNA_protein= (mass_dry * 0.75)

# mass_RNA / mass_protein = ( mass_RNA_protein * frac_RNA ) / ( mass_RNA_protein * frac_protein )
# * frac_RNA) / (mass_dry * 0.75 * frac_protein)

mass_RNA = RNA_protein_ratio * (mass_dry * 0.75) / (1+ RNA_protein_ratio)
mass_protein =(mass_dry * 0.75)/ (1+RNA_protein_ratio)



fig, ax = plt.subplots(1, 2, figsize=(8,4))

ax[0].plot(growth_rate, mass_dry, label = 'dry mass')
ax[0].plot(growth_rate, mass_protein, label = 'protein')
ax[0].plot(growth_rate, mass_RNA, label = 'RNA')
ax[0].set_xlabel(r'growth rate $\lambda$')
ax[0].set_ylabel('mass (fg)')
ax[0].legend()

ax[1].plot(growth_rate, mass_dry/vol, label = 'mass conc.')
ax[1].plot(growth_rate, mass_protein/vol, label = 'protein conc.')
ax[1].plot(growth_rate, mass_RNA/vol, label = 'RNA conc.')
ax[1].set_xlabel(r'growth rate $\lambda$')
ax[1].set_ylabel('concentration (fg/fL)')

ax[1].set_ylim(0,340)

fig.savefig('../../../figures/mass_predictions.pdf')
