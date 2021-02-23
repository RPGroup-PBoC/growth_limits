#%%
import pandas as pd 
import numpy as np 

# # Load the compiled data
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')

L_R = 7459.0 # length of all subunits in ribosomes, in amino acids

ribosome_genes = ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE',
              'rpsF', 'rpsG', 'rpsH', 'rpsI', 'rpsJ', 'rpsK',
              'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
              'rpsR', 'rpsS', 'rpsT', 'rpsU', 'sra', 'rplA', 'rplB',
              'rplC', 'rplD', 'rplE', 'rplF', 'rplJ',
              'rplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP', 'rplQ',
              'rplR', 'rplS','rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY',
              'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', 'rpmF', 'rpmG', 'rpmH',
              'rpmI', 'rpmJ', 'ykgM', 'ykgO']

# Add in data from Dai et al.
dai_df = pd.DataFrame(np.array(\
        [[0.0, 0.08853279947047754],
        [0.03324706890845652, 0.09742834356027935],
        [0.12844176066233703, 0.12157153165470358],
        [0.19652012565308674, 0.12917148174069676],
        [0.23055930814846148, 0.13297145678369332],
        [0.2849547247405284, 0.14694954887503597],
        [0.33601892391911736, 0.15201256530868013],
        [0.4074068045812377, 0.17107537557577435],
        [0.417639175984852, 0.16979497279144085],
        [0.4517000602223341, 0.17104716331103476],
        [0.485674137491387, 0.18249049192423916],
        [0.5503561798423366, 0.1888187199227418],
        [0.6727865579409387, 0.21549233114688282],
        [0.6864152519843529, 0.21548365045003987],
        [0.7000547968988209, 0.21420107749149564],
        [0.7170798135820351, 0.21546411888214323],
        [0.744196140345166, 0.2320073568905744],
        [0.9177883754618401, 0.25227895419304786],
        [0.9448830004828637, 0.27136997672488156],
        [0.9926268331190284, 0.2662440252391261],
        [0.9753630972726335, 0.29300661360590724],
        [1.0979236858238794, 0.3043935176896325],
        [1.1624538159800777, 0.32855623735195344],
        [1.2677832212980895, 0.36288405301735593],
        [1.566952587118931, 0.4404005056505911],
        [1.7949076862145108, 0.4784718718295111]]),
        columns=['growth_rate_hr', 'mass_fraction'])
dai_df['source'] = 'Dai et al., 2016'
# add in data from Scott et al.
scott_df = pd.DataFrame(np.array(\
        [[0.4, 0.177],
        [0.57, 0.230],
        [0.71, 0.221],
        [1.00, 0.287],
        [1.31, 0.414],
        [1.58, 0.466]]), columns=['growth_rate_hr', 'mass_fraction'])
scott_df['source'] = 'Scott et al., 2010'
# Add in data from Forchhammer & Lindahl
forchammer_df = pd.DataFrame(np.array(\
        [[0.38, 0.189],
        [0.60, 0.224],
        [1.04, 0.295],
        [1.46, 0.421],
        [1.73, 0.469]]), columns=['growth_rate_hr', 'mass_fraction'])
forchammer_df['source'] = 'Forchammer & Lindahl, 1971'

# Bremmer and Dennis
bremmer_df = pd.DataFrame(np.array(\
        [[0.42, 0.200],
        [0.69, 0.225],
        [1.04, 0.331],
        [1.39, 0.391],
        [1.73, 0.471]]),
        columns=['growth_rate_hr', 'mass_fraction'])
bremmer_df['source'] = 'Bremmer & Dennis, 2008'

df_ribo_frac = pd.DataFrame()
for c, d in data.groupby(['dataset', 'condition', 'growth_rate_hr']):
    mass_ribo = d[d['gene_name'].isin(ribosome_genes)].fg_per_cell.sum()
    frac_ribo = (mass_ribo )/ d.fg_per_cell.sum()

    data_list = {'mass_fraction' : frac_ribo,
                'dataset' : c[0],
                'condition' : c[1],
                'growth_rate_hr' : c[2]}

    df_ribo_frac = df_ribo_frac.append(data_list,
                                        ignore_index = True)
df_ribo_frac = df_ribo_frac[(df_ribo_frac['dataset'] != 'valgepea_2013') &
                            (df_ribo_frac['dataset'] != 'peebo_2015')]
df_ribo_frac.loc[df_ribo_frac['dataset']=='li_2014', 'source'] = 'Li et al., 2014'
df_ribo_frac.loc[df_ribo_frac['dataset']=='schmidt_2016', 'source'] = 'Schmidt et al., 2016'
df_ribo_frac = df_ribo_frac[['growth_rate_hr', 'mass_fraction', 'source']]

# Concatenate and save everything
df = pd.concat([dai_df, bremmer_df, forchammer_df, scott_df], sort=False)
df['mass_fraction'] = 0.4558 * df['mass_fraction'].values
df = pd.concat([df_ribo_frac, df], sort=False)
df.to_csv('../../data/mass_fraction_compiled.csv', index=False)
df
# %%
