#%%
import pandas as pd 
import tqdm
import numpy as np
import prot.size

# Load  the necessary datasets. 
data = pd.read_csv('../../data/compiled_annotated_complexes.csv', comment='#')
data.dropna(subset=['n_units'], inplace=True)

# define necessary complexes. 
complexes = {'dnap': {'name': 'DNA polymerase III (holo enzyme)', 
                       'complexes': ['CPLX0-3803'],
                       'rate_per_sec': 600,
                       'units': 'bp/s',
                       'method': 'sum',
                       'category':'synthesis'},
             'rnap': {'name':'RNA polymerase (core enzyme)',
                       'complexes': ['APORNAP-CPLX'],
                       'rate_per_sec': 40,
                       'units':'nt/s',
                       'method':'sum',
                       'category':'synthesis'},
            'dntp': {'name': 'Ribonucleoside-diphosphate reductase (I)',
                     'complexes': ['RIBONUCLEOSIDE-DIP-REDUCTI-CPLX',
                                   'RIBONUCLEOSIDE-DIP-REDUCTII-CPLX',
                                   'NRDACTMULTI-CPLX'],
                     'rate_per_sec': 10,
                     'units': 'dNTP/s' ,
                     'method': 'sum',
                     'category':'synthesis'},
            'sigma70': {'name':'σ-70 (RpoD)',
                       'gene_name': ['rpoD'],
                       'rate_per_sec': 40,
                       'units':'nt/s',
                       'method':'sum',
                       'category':'synthesis'},
            'all_sigma': {'name': 'all σ-factors',
                          'gene_name': ['rpoE', 'fecI', 'rpoF', 'rpoH', 'rpoN', 'rpoD', 'rpoS'],
                          'rate_per_sec': 40,
                          'units':'nt/s',
                          'method':'sum',
                          'category':'synthesis'},
             'trna': {'name': 'tRNA ligases', 
                      'gene_name': ['argS', 'cysS', 'glnS', 'gltX', 
                                  'ileS', 'leuS', 'valS', 'alaS',
                                  'asnS', 'aspS', 'tyrS', 'trpS',
                                  'thrS', 'serS', 'proS', 'pheS',
                                  'metG', 'lysS', 'hisS', 'glyS'],
                    'rate_per_sec': 20,
                    'units': 'AA/s',
                    'method':'sum',
                    'category':'synthesis'},
            'carbon_tport': {'name': 'Carbohydrate Transporters (average)', 
                              'go_terms': ['GO:0009401'],
                              'rate_per_sec': 200,
                              'units': 'carbs/s',
                              'method': 'avg',
                              'category':'transport'},
            'carbon_tport_tot': {'name': 'Carbohydrate Transporters (total)', 
                              'go_terms': ['GO:0009401'],
                              'rate_per_sec': 200,
                              'units': 'carbs/s',
                              'method': 'sum',
                              'category':'transport'},
            'glucose_tport': {'name': 'Glucose/Mannose Transporters',
                              'complexes': ['CPLX-157', 'CPLX-165'],
                              'rate_per_sec': 200,
                              'units': 'gluc/sec',
                              'method':'sum',
                              'category':'transport'},
             'glycerol_tport': {'name': 'Glycerol Transporters',
                              'complexes': ['CPLX0-7654'],
                              'rate_per_sec': 0,
                              'units': 'glyc/sec',
                              'method':'sum',
                              'category':'transport'},
            'galactose_tport': {'name':'Galactose Transporter (MglABC)',
                                'complexes': ['ABC-18-CPLX'],
                                'rate_per_sec': 10,
                                'units': 'gal/s',
                                'method':'sum',
                                'category':'transport'},
            'xylose_tport': {'name': 'Xylose Transporters',
                              'gene_name': ['xylG', 'xylH', 'xylF', 'xylE'],
                              'rate_per_sec': 200,
                              'units': 'xyl/sec',
                              'method': 'sum',
                              'category':'transport'},
            'glucosamine_tport': {'name': 'Glucosamine Transporter NagE',
                                  'complexes': ['CPLX-167'],
                                  'rate_per_sec': 200,
                                  'units': 'glcn/s',
                                  'method':'sum',
                                  'category':'transport'},
            'fructose_tport': {'name': 'Fructose Transporter FruBA',
                                  'complexes': ['CPLX-158'],
                                  'rate_per_sec': 200,
                                  'units': 'frc/s',
                                  'method':'sum',
                                  'category':'transport'},
            'nitrogen_tport': {'name': 'Ammonium Transporter (AmtB)',
                               'gene_name': ['amtB'],
                               'rate_per_sec': 300,
                               'units': 'NH4+/s',
                               'method':'sum',
                               'category':'transport'},
            'sulfur_tport': {'name': 'Sulfate Transporter (CysUWA)',
                            'complexes': ['ABC-70-CPLX', 'ABC-7-CPLX'],
                            'rate_per_sec':10,
                            'units': 'SO4/s',
                            'method':'avg',
                            'category':'transport'},
            'phosphate_tport': {'name': 'Phosphate Transport System',
                            'gene_name': ['pitA', 'pitB'],
                            'rate_per_sec': 0.01,
                            'units': 'Pi/s',
                            'method':'sum',
                            'category':'transport'},
            'ribosome': {'name': 'Ribosome (50S + 30S)',
                         'complexes': ['CPLX0-3964'],
                         'rate_per_sec': 15,
                         'units': 'AA/s',
                         'method':'sum',
                         'category':'synthesis'},
            'eftu': {'name': 'Elongation Factor EF-Tu', 
                     'gene_name': ['tufA', 'tufB'],
                     'method': 'sum',
                     'rate_per_sec': 20,
                     'units':'peptide bonds/s',
                     'category':'synthesis'},
            'atp_synthase': {'name': 'F1-F0 ATP Synthase',
                            'complexes':  ['ATPSYN-CPLX'],
                            'method': 'sum', 
                            'rate_per_sec': 300,
                            'units':'atp/s',
                            'category': 'energy production'}, 
            'proton_gradient': {'name': 'respiratory complex', 
                    'go_terms': ['GO:0019646', 'GO:0006136', 'GO:0006137', 'GO:0006138'],
                    'method':'sum',
                    'rate_per_sec': 5E3,
                    'units': 'protons/s',
                    'category': 'energy production'},
            'fas': {'name': 'Fatty Acid Synthetases (FabA + FabZ)',
                    'gene_name': ['fabZ', 'fabA'],
                    'method': 'sum', 
                    'rate_per_sec': 1,
                    'units':'lipid/s',
                    'category': 'synthesis'}}

# %%
complex_df = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr'])):
    for k, v in complexes.items():
        if 'complexes' in list(v.keys()):
            _d = d[d['complex'].isin(v['complexes'])]
        if 'go_terms' in list(v.keys()):
            cplxs = []
            for t in v['go_terms']:
                __d = d[d['go_terms'].str.contains(t)]
                for kplx in __d['complex'].unique():
                    cplxs.append(kplx)
            _d = d[d['complex'].isin(cplxs)]
        if 'gene_name' in list(v.keys()):
            _d = d[d['gene_name'].isin(v['gene_name'])]
        if len(_d) > 0:
            _d = _d.drop_duplicates(subset=['gene_name'])
            _d = _d.groupby(['complex'])['n_units'].mean().reset_index()
            if v['method'] == 'sum':
                units = _d['n_units'].sum()
                _method = 'sum total'
            if v['method'] == 'avg':
                units = _d['n_units'].mean()
                _method = 'average'

            # assemble a dictionary 
            volume = np.round(prot.size.lambda2size(g[3]), 2)
            _data = {'dataset':g[0], 'dataset_name':g[1],
                     'condition':g[2], 'growth_rate_hr':g[3],
                     'volume': volume,
                     'n_complex':units, 
                     'rate': v['rate_per_sec'],
                     'rate_units': v['units'],
                     'shorthand': k,
                     'name': v['name'],
                     'aggregation_method': _method,
                     'category': v['category'],
                     'concentration_uM': 1E6 * (units / 6.022E23) / (volume * 1E-15)}
        
            complex_df = complex_df.append(_data, ignore_index=True)

complex_df.to_csv('../../data/compiled_estimate_categories.csv', index=False)
# %%







# %%
