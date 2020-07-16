#%%
import numpy as np
import pandas as pd

##############################################
##############################################
# gather cell size information and save to csv
##############################################
##############################################
size_df = pd.DataFrame()

# Taheri-Araghi et al. 2015; X = # divisions per hr (so need to divide by ln(2) to get lambda)
                               # Y = cell size/ volume in um^3

[division_rate_TA, vol_TA] = np.array(\
            [[1.16, 0.612],
            [1.17, 0.653],
            [1.57, 0.866],
            [1.97, 1.16],
            [2.23, 1.50],
            [2.64, 2.38],
            [3.48, 3.94]]).T

lambda_TA = division_rate_TA * np.log(2)

size_df = size_df.append(pd.DataFrame({'volume_um3' : vol_TA,
                            'growth_rate_hr' : lambda_TA,
                            'dataset_name' : 'Taheri-Araghi et al 2015'},
                    columns = ['volume_um3', 'growth_rate_hr', 'dataset_name']))

# Basan et al. 2015   X = growth rate lambda (hr-1), Y = cell size/volume (um^3)
[lambda_basan, vol_basan] = np.array(\
             [[0.42, 1.37],
            [0.47, 1.70],
            [0.70, 1.70],
            [0.98, 2.34],
            [1.29, 3.34],
            [1.84, 6.62]]).T

size_df = size_df.append(pd.DataFrame({'volume_um3' : vol_basan,
                            'growth_rate_hr' : lambda_basan,
                            'dataset_name' : 'Basan et al 2015'},
                    columns = ['volume_um3', 'growth_rate_hr', 'dataset_name']))

# load Schmidt 2016 data
schmidt_sizes_df = pd.read_csv('../../../data/schmidt2016_raw_data/schmidt2016_growth_rates.csv')

size_df = size_df.append(pd.DataFrame({'volume_um3' : schmidt_sizes_df.volume_fL.values,
                            'growth_rate_hr' : schmidt_sizes_df.growth_rate_hr.values,
                            'dataset_name' : 'Schmidt et al 2016'},
                    columns = ['volume_um3', 'growth_rate_hr', 'dataset_name']))


# Si et al. 2017, X = growth rate lambda (hr-1), Y = cell size/volume (um^3)
[lambda_Si, vol_Si] = np.array(\
           [[0.36, 0.45],
            [0.62, 0.61],
            [0.64, 0.61],
            [0.64, 0.60],
            [0.70, 0.57],
            [0.76, 0.53],
            [0.77, 0.51],
            [1.05, 0.89],
            [1.10, 0.98],
            [1.03, 1.04],
            [1.02, 1.043],
            [1.10, 1.09],
            [1.11, 1.29],
            [1.16, 1.30],
            [1.16, 1.24],
            [1.26, 1.37],
            [1.32, 1.35],
            [1.35, 1.18],
            [1.51, 1.42],
            [1.51, 1.62],
            [1.51, 1.62],
            [1.56, 1.62],
            [1.81, 2.75],
            [1.77, 3.26],
            [1.71, 3.26],
            [1.71, 3.42],
            [1.76, 3.90],
            [1.94, 3.26],
            [1.95, 3.47],
            [2.12, 4.18]]).T

size_df = size_df.append(pd.DataFrame({'volume_um3' : vol_Si,
                            'growth_rate_hr' : lambda_Si,
                            'dataset_name' : 'Si et al 2017'},
                    columns = ['volume_um3', 'growth_rate_hr', 'dataset_name']))



size_df.to_csv('../../../data/cell_size_literature/size_data.csv')
