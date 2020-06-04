#%%
import numpy as np
import pandas as pd
from scipy import stats
import tqdm
import prot.size as size

# Load the data quantifying absolute protein synthesis rates.
counts = pd.read_csv('../../../data/peebo2015_raw_data/peebo2014_copynums_minimal.csv')
counts = counts.append(pd.read_csv('../../../data/peebo2015_raw_data/peebo2014_copynums_rich.csv'))

# Get a list of the unique conditions.
conditions = counts['condition'].unique()[:2]

# Load the annotation list.
colicogs = pd.read_csv('../../../data/ecoli_genelist_master.csv')

# Instantiate a blank dataframe to which measurements and annotations will be
# added
dfs = []

# Iterate through each gene in the count data.
for g, d in tqdm.tqdm(counts.groupby('gene'), desc="Iterating through genes..."):

    # Determine number of entries per gene
    gene = colicogs[colicogs['gene_name'].str.lower()==g.lower()]
    b_number = gene['b_number'].unique()[0]
    gene_product = gene['gene_product'].unique()[0]
    go_term = ';'.join(list(gene['go_terms'].unique()))
    if len(gene) > 0:
        cog_class = gene['cog_class'].values[0]
        cog_cat = gene['cog_category'].values[0]
        cog_letter = gene['cog_letter'].values[0]
        gene_product= gene['gene_product'].values[0]
        mw = gene['mw_fg'].values[0]
        go_term = ';'.join(list(gene['go_terms'].unique()))
        for _c, _d in d.groupby(['growth_rate_hr-1']):
            # volume predictions based on MG1655 data, Si, F. et al. (2017, 2019)
            vol = size.lambda2size(_c)
            # extract relevant information.
            gene_dict = {
                'gene_name': g.lower(),
                'b_number': b_number,
                'condition': _d['condition'].unique()[0],
                'corrected_volume': vol,
                'reported_tot_per_cell': _d['copy_number_molecule-per-fL'].values[0] * vol,
                'reported_fg_per_cell': _d['copy_number_molecule-per-fL'].values[0] * vol * mw,
                'go_terms':go_term,
                'cog_class': cog_class,
                'cog_category': cog_cat,
                'cog_letter': cog_letter,
                'gene_product': gene_product,
                'growth_rate_hr': _c
                }
            dfs.append(pd.DataFrame(gene_dict, index=[0]))
    else:
        print(f'Warning!!! {g} not found in the gene list! Not including in final tally')

df = pd.concat(dfs, sort=False)

#%%
# Iterate through the conditions and correct fg and tot_per_cell as necessary.
# for g in df['growth_rate_hr'].unique():
#     gr = g
#     rel_corr_fg = df[df['growth_rate_hr'] == g]['reported_fg_per_cell'].sum() / \
#             size.lambda2P(gr)
#
#     print(g, ': total mass fg: ', size.lambda2P(gr),
#             ' volume: ', size.lambda2size(gr),
#             ' relative change in total fg: ', 1/rel_corr_fg)
#
#     df.loc[df['growth_rate_hr']==g, 'tot_per_cell'] = \
#                     df.loc[df['growth_rate_hr']==g]['reported_tot_per_cell'] / rel_corr_fg
#     df.loc[df['growth_rate_hr']==g, 'fg_per_cell'] =  \
#                     df.loc[df['growth_rate_hr']==g]['reported_fg_per_cell'] / rel_corr_fg

for g, d in df.groupby(['condition', 'growth_rate_hr']):
    gr = g[1]
    rel_corr_fg = df[(df['condition'] == g[0]) & \
                    (df['growth_rate_hr'] == g[1])]['reported_fg_per_cell'].sum() / \
            size.lambda2P(gr)

    print(g, ': total mass fg: ', np.round(size.lambda2P(gr),2),
            ' volume: ', np.round(size.lambda2size(gr),2),
            ' relative change in total fg: ', np.round(1/rel_corr_fg, 2))

    df.loc[(df['condition']==g[0]) &
            (df['growth_rate_hr']==g[1]), 'tot_per_cell'] = df.loc[(df['condition']==g[0]) &
                    (df['growth_rate_hr']==g[1])]['reported_tot_per_cell'] / rel_corr_fg
    df.loc[(df['condition']==g[0]) &
            (df['growth_rate_hr']==g[1]), 'fg_per_cell'] =  df.loc[(df['condition']==g[0]) &
                    (df['growth_rate_hr']==g[1])]['reported_fg_per_cell'] / rel_corr_fg


#%%
# include dataset notation.
df['dataset'] = 'peebo_2015'
df['dataset_name'] = 'Peebo et al. 2015'
df['strain'] = 'BW25113'

#%%
df.to_csv('../../../data/peebo2015_longform_annotated.csv')
# %%

##########################################
# #%% Ignore section for now; use correction above.
#%%
# The reported total fg/fL is lower than reported in Schmidt and Valgepea;
# Schmidt and Valgepea appear to have roughly consistent fg/fL;
# Due to the the discrepency, and the lack of cell counts to calculate
# copies or fg -- per cell --, we are imposing two constraints in order
# to calculate per cell quantitites.
# 1. The total fg at each growth rate should be consistent with Schmidt
# (which we believe to be the most comprehensive in their quantification of total
# mass per cell). Perform regression on Schmidt fg vs. growth rate and use this
# to predict the mass per cell for each Peebo condition.
# 2. (DONE ABOVE) Calculate fg per cell using estimates of cell volume from
# the 'growth law' relations quantified in Si, F. et al. (2017), Current Biology.
# 3. Apply correction factor to Peebo quantities; this'll effectively adjust
# their values so that the fg/fL is consistent with Schmidt and Valgepea.

# Determine fg as a function of growth rate using Schmidt et al. data (which
# we have also corrected for discrepencies in cell volume, now using the
# volume predictions of Si, F. et al. (2017), Current Biology.

# # #%% STEP 1: Determine expected fg as function of growth rate using Schimdt data
# bnum_list = df.b_number.unique()
# df_schmidt = pd.read_csv('../../../data/schmidt2016_longform_annotated.csv')
# df_schmidt = df_schmidt[df_schmidt['b_number'].isin(bnum_list)]
# df_schmidt = df_schmidt[['b_number',  'condition', 'growth_rate_hr',
#                             'corrected_volume', 'fg_per_cell']]
# df_schmidt = df_schmidt[df_schmidt.growth_rate_hr > 0]
# df_schmidt = df_schmidt.drop_duplicates().sort_values(by='growth_rate_hr')
#
# # for cond, data in df_schmidt.groupby(['condition']):
# #     fg_perfL = data.fg_per_cell.sum()/ data.corrected_volume.unique()[0]
# fg = [data.fg_per_cell.sum() for
#                     cond, data in df_schmidt.groupby(['condition'], sort=False)]
# growth_rate_hr = [data.growth_rate_hr.unique()[0] for
#                     cond, data in df_schmidt.groupby(['condition'], sort=False)]
#
# # Fit line for fg/fL as function of growth rate
# # i.e. create a linear regression model
# slope, intercept, r_value, p_value, std_err = stats.linregress(growth_rate_hr, fg)
# # report coefficient of determination,  r**2
# print("r-squared:", r_value**2)
#
# # #%% STEP 3: Apply correction to Peebo dataset
#
# # # Compute the mass corrections
# _conditions = df[['gene_name', 'condition', 'reported_fg_per_cell', 'growth_rate_hr']].drop_duplicates()
# _conditions = _conditions.groupby(['growth_rate_hr', 'condition']).sum().reset_index()
# _conditions['correction_factor'] = _conditions['reported_fg_per_cell'].values / (intercept + slope * _conditions['growth_rate_hr'].values )
#
# # Apply corrections to counts/ fg per cell
# for g, d in _conditions.groupby(['growth_rate_hr', 'condition']):
#     corr_factor = _conditions[(_conditions['growth_rate_hr']==g[0]) & (_conditions['condition']==g[1])]['correction_factor'].values[0]
#     df.loc[(df['growth_rate_hr']==g[0]) & (df['condition']==g[1]), 'tot_per_cell'] = \
#             df.loc[(df['growth_rate_hr']==g[0]) & (df['condition']==g[1])]['reported_tot_per_cell'] / corr_factor
#     df.loc[(df['growth_rate_hr']==g[0]) & (df['condition']==g[1]), 'fg_per_cell'] =  \
#             df.loc[(df['growth_rate_hr']==g[0]) & (df['condition']==g[1])]['reported_fg_per_cell'] / corr_factor
##########################################




# %%
