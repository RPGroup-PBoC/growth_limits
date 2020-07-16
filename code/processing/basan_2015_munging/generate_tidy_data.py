#%%
import numpy as np
import pandas as pd

##############################################
##############################################
# gather DNA mass data in Basan 2015; from SI PDF
##############################################
##############################################

# growth rate (hr-1), dna (femtogram per cell)
lambda_basan_dna = [0.42, 0.45, 0.7, 0.98, 1.27, 1.84]
dna_basan = 1E15*np.divide(1E-6*np.array([16.5, 14.2, 14.1, 11.9, 11.3, 11.1]),
                      1E8*np.array([19.4, 17.1, 16.0, 10.7, 7.93, 3.43]))

#  protein data from appendix Figure S7 for nutrient limitation
protein_basan = [202,  223, 224, 315, 405, 965]


dai_df = pd.DataFrame({
     'growth_rate_hr' : lambda_basan_dna,
     'dna_fg' : dna_basan,
     'protein_fg' : protein_basan},
        columns = ['growth_rate_hr',
                   'dna_fg', 'protein_fg'])

dai_df.to_csv('../../../data/basan2015_raw_data/basan2015_data.csv')
