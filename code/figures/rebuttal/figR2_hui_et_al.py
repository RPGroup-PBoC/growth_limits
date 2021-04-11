#%%
import numpy as np
import pandas as pd
import glob as glob

import matplotlib.pyplot as plt
from collections import OrderedDict
import prot.viz
import prot.estimate
colors = prot.viz.plotting_style()



data = pd.read_csv('../../../data/compiled_absolute_measurements.csv')
#compare with Schmidt glucose
data = data[data.dataset == 'schmidt_2016']
data = data[data.condition != '42C']


# import panel as pn
# import altair as alt
#

# import tqdm
# import prot.viz
# prot.viz.altair_theme()
# pn.extension('vega')
# _ = alt.data_transformers.enable('default')
# alt.data_transformers.disable_max_rows()

gr_dict = {'NQ381_Lac_3MBA0': np.log(2)/(92/60), 'NQ381_Lac_3MBA25': np.log(2)/(72/60),
             'NQ381_Lac_3MBA50': np.log(2)/(62/60), 'NQ381_Lac_3MBA500': np.log(2)/(48/60),
           'NCM3722_Lac_3MBA0': np.log(2)/(40/60),
             'NCM3722_Glc_Cm8': np.log(2)/(147/60), 'NCM3722_Glc_Cm4': np.log(2)/(102/60),
             'NCM3722_Glc_Cm2': np.log(2)/(65/60), 'NCM3722_Glc_Cm0': np.log(2)/(42/60),
          'NQ393_Glc_IPTG30': np.log(2)/(91/60),  'NQ393_Glc_IPTG40': np.log(2)/(69/60),
                                         'NQ393_Glc_IPTG50': np.log(2)/(58/60),
          'NQ393_Glc_IPTG100': np.log(2)/(47/60), 'NCM3722_Glc_IPTG0': np.log(2)/(43/60)}





# df_ = pd.read_csv('../../data/hui_2015_raw/header_removed/hui_2015_merged_proteomic.csv')
dir_ = '../../../data/hui_2015_raw/header_removed/'
a_sector = pd.read_csv(dir_ + 'Hui_2015_a_sector_genes.csv').gene_name.values
r_sector = pd.read_csv(dir_ + 'Hui_2015_r_sector_genes.csv').gene_name.values
o_sector = pd.read_csv(dir_ + 'Hui_2015_o_sector_genes.csv').gene_name.values
s_sector = pd.read_csv(dir_ + 'Hui_2015_s_sector_genes.csv').gene_name.values
c_sector = pd.read_csv(dir_ + 'Hui_2015_c_sector_genes.csv').gene_name.values
u_sector = pd.read_csv(dir_ + 'Hui_2015_u_sector_genes.csv').gene_name.values

df_frac = pd.read_csv('../../../data/hui_2015_raw/header_removed/hui_2015_merged_proteomic_fractional_abundances.csv')

df_sector_frac = pd.DataFrame()
for c, d in df_frac.groupby(['exp_type', 'experiment', 'trial']):

#     frac_ribo = d[d['gene'].isin(ribosome_genes)].frac_gene.sum()
    frac_a_sector = d[d['gene'].isin(a_sector)].frac_gene.sum()
    frac_c_sector = d[d['gene'].isin(c_sector)].frac_gene.sum()
    frac_o_sector = d[d['gene'].isin(o_sector)].frac_gene.sum()
    frac_r_sector = d[d['gene'].isin(r_sector)].frac_gene.sum()
    frac_s_sector = d[d['gene'].isin(s_sector)].frac_gene.sum()
    frac_u_sector = d[d['gene'].isin(u_sector)].frac_gene.sum()

    data_list = {'frac_a_sector' : frac_a_sector,
                 'frac_c_sector' : frac_c_sector,
                 'frac_o_sector' : frac_o_sector,
                 'frac_r_sector' : frac_r_sector,
                 'frac_s_sector' : frac_s_sector,
                 'frac_u_sector' : frac_u_sector,
                 'exp_type' : c[0],
                'experiment' : c[1],
                'trial' : c[2]}

    df_sector_frac = df_sector_frac.append(data_list,
                                        ignore_index = True)


source_data = pd.DataFrame()

for c, d in df_sector_frac.groupby(['exp_type', 'experiment', 'trial']):
#     if '' in c[1]:
#         continue
#     frac_ribo = d[d['gene'].isin(ribosome_genes)].frac_gene.sum()
    frac_a_sector_avg = d.frac_a_sector.mean()

    frac_c_sector_avg = d.frac_c_sector.mean()
    frac_o_sector_avg = d.frac_o_sector.mean()
    frac_r_sector_avg = d.frac_r_sector.mean()
    frac_s_sector_avg = d.frac_s_sector.mean()
    frac_u_sector_avg = d.frac_u_sector.mean()

    data_list = {'frac_a_sector' : frac_a_sector_avg,
                 'frac_c_sector' : frac_c_sector_avg,
                 'frac_o_sector' : frac_o_sector_avg,
                 'frac_r_sector' : frac_r_sector_avg,
                 'frac_s_sector' : frac_s_sector_avg,
                 'frac_u_sector' : frac_u_sector_avg,
                 'exp_type' : c[0],
                'experiment' : c[1],
                'growth_rate_hr' : gr_dict[c[1]]}

    source_data = source_data.append(data_list,
                                        ignore_index = True)

source_data_ = source_data[source_data.exp_type == 'C-lim']
source_data_ = source_data_.append(source_data[source_data.exp_type == 'A-lim'])
source_data_ = source_data_.append(source_data[source_data.exp_type == 'R-lim'])


fig, ax = plt.subplots(1, 2, figsize = (8,4))
ax = ax.ravel()


######################

lim_colors = {'C-lim' : colors['red'], 'A-lim' : colors['blue'], 'R-lim' : colors['green']}

for c, d in source_data_.groupby(['exp_type', 'experiment']):
    ax[0].errorbar(d.growth_rate_hr.unique(), d.frac_r_sector.mean(),
                   yerr=d.frac_r_sector.std(), color = lim_colors[c[0]],
                  fmt = 'o', label = c[0])
ax[0].legend(fontsize = 10)
ax[0].set_ylim(0,0.5)
ax[0].set_xlim(0,1.3)
ax[0].set_title('R-sector', fontsize = 12)
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize = 12)
ax[0].set_ylabel('proteomic fraction', fontsize = 12)

handles, labels = ax[0].get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax[0].legend(by_label.values(), by_label.keys(), loc = 'upper right', fontsize = 10)

#####################

for i, cond in enumerate(['NCM3722_Glc_IPTG0']):
    # print(cond)
    # if cond != 'NCM3722_Glc_IPTG0':
    #     continue

    df_hui = df_frac[df_frac.experiment == cond]
    df_schmidt = data[data.condition == 'glycerol_pAA']
    df_schmidt['rel'] = df_schmidt.tot_per_cell/df_schmidt.tot_per_cell.sum()

    df_compare = pd.DataFrame()
    for gene, d in df_schmidt.groupby('gene_name'):
        if gene in df_hui.gene.values:
            data_list = {'schmidt_abs' : d.tot_per_cell.values[0],
                         'schmidt' : d.rel.values[0],
                     'hui' : df_hui[df_hui.gene == gene].frac_gene.mean()}
            df_compare = df_compare.append(data_list,
                                          ignore_index = True)

    #  plot on log scale; ignore zero values
    df_compare = df_compare[df_compare.hui != 0]
    df_compare = df_compare[df_compare.schmidt != 0]


    ax[1].plot(df_compare['schmidt'].values, df_compare['hui'].values, 'o', ms=2,
            alpha=0.5, markeredgewidth=0, markeredgecolor='k')
    ax[1].set_xlabel('relative protein abundance\n(Schmidt,' + 'glycerol_pAA' +  ')', fontsize = 12)
    ax[1].set_ylabel('relative protein abundance\n(Hui,' + cond +  ')', fontsize = 12)

    ax2 = ax[1].twiny()
    ax2.plot(df_compare['schmidt_abs'].values, df_compare['hui'].values, 'o', ms=2,
            alpha=0.5, markeredgewidth=0, markeredgecolor='k')
    ax2.lines = []
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    #I liked the log-scaled x-ticks that got generated with the twiny
    # lets force this for the y axis by  repeating
    ax2_ = ax[1].twinx()
    ax2_.set_yticks([])

for _ax in ax[1:]:
    _ax.set_xscale('log')
    _ax.set_yscale('log')
    _ax.set_xlim(1E-6,1E-1)
    _ax.set_ylim(1E-5,1E-1)
    _ax.tick_params(axis='x', labelsize=10)
    _ax.tick_params(axis='y', labelsize=10)


plt.tight_layout()

plt.savefig('../../../figures/figR2_hui_data_comparison.pdf')
