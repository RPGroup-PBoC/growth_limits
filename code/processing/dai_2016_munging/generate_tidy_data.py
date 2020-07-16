#%%
import numpy as np
import pandas as pd

##############################################
##############################################
# gather grab nutrient limitation data from Dai et al. 2016
#  data from SI PDF
##############################################
##############################################

conditions = ['RDM + 0.2% glucose+10 mM NH4Cl',
            '0. 2 % glucose+cAA+10 mM NH4Cl',
           ' 10 mM glucose-6-phosphate+10 mM gluconate +10 mM NH4Cl',
            '0.2% glucose+10 mM NH4Cl',
            '0.2% xylose+10 mM NH4Cl',
            '0.2 % glycerol+10 mM NH4Cl',
            '0.2% fructose+10 mM NH4Cl',
            '0.2% sorbitol+10 mM NH4Cl',
            '0.2% galactose+10 mM NH4Cl',
            '60 mM acetate+10 mM NH4Cl',
            '0.2% mannose+10 mM NH4Cl',
            '0.1% mannose+10 mM NH4Cl',
            '20 mM potassium aspartate',
            '0.075% mannose+10 mM NH4Cl',
            '20 mM aspartate+10 mM NH4Cl',
            '0.2% glycerol +10 mM Arginine',
            '20 mM glutamate+10 mM NH4Cl',
            '0.2% glycerol+20 mM Threonine']

# measured elongation rate using LacZ induction assay
erate = [16.7, 16.3, 16.1, 15.9, 14.9, 15.0, 14.7, 13.7, 13.1, 12.6, 13.0, 12.4,
        12.0, 12.1, 12.3, 11.6, 10.7, 9.4]

# total RNA-to-protein ratio
RNA_P = [0.476, 0.364, 0.306, 0.294, 0.233, 0.227, 0.217, 0.193, 0.184, 0.172,
        0.172, 0.152, 0.152, 0.147, 0.137, 0.130, 0.118, 0.097]

# ribosomal fraction from mass spec. experiments (R-protein fraction based on sum of all
# ribosomal subunits)
r_prot_frac = [np.nan, np.nan, np.nan, 11.6, np.nan, np.nan, np.nan, np.nan,
        np.nan, 7.2, np.nan, np.nan, np.nan, 5.1, np.nan, 6.1, 4.7, 4.4]

# growth rate, fraction active ribosomes
gr, fa = np.array([[1.8, 0.958],
            [1.28, 0.9],
            [1.12, 0.927],
            [0.98, 0.865],
            [0.75, 0.902],
            [0.69, 0.849],
            [0.69, 0.888],
            [0.55, 0.879],
            [0.5, 0.860],
            [0.46, 0.879],
            [0.41, 0.756],
            [0.34, 0.751],
            [0.33, 0.756],
            [0.29, 0.683],
            [0.23, 0.590],
            [0.201, 0.554],
            [0.13, 0.441],
            [0.035, 0.168]]).T

dai_df = pd.DataFrame({'condition' : conditions,
     'Cm (μM)' : np.zeros(len(gr)),
     'RNA_P_ratio' : RNA_P,
     'growth_rate_hr' : gr,
     'Translational elongation rate (aa/s)' : erate,
     'measured_prot_frac' : r_prot_frac,
     'f_a' : fa,
    'type' : ['nutrient limitation' for i in np.arange(len(gr))]},
        columns = ['condition', 'Cm (μM)', 'RNA_P_ratio', 'growth_rate_hr',
                   'Translational elongation rate (aa/s)', 'measured_prot_frac',
                   'f_a', 'type'])

dai_df.to_csv('../../../data/dai2016_raw_data/dai2016_summary.csv')
