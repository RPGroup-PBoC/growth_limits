#%%
import pandas as pd 
import tqdm

# Load  the necessary datasets. 
data = pd.read_csv('../../data/compiled_annotated_complexes.csv', comment='#')

# define necessary complexes. 
complexes = {'dnap': {'name': 'DNA polymerase III (core enzyme)', 
                       'complexes': ['CPLX0-3803'],
                       'rate_per_sec': 600,
                       'units': 'bp/s',
                       'method': 'sum'},
             'rnap': {'name':'RNA polymerase (core enzyme)',
                       'complexes': ['APORNAP-CPLX'],
                       'rate_per_sec': 20,
                       'units':'nt/s',
                       'method':'sum'},
             'trna': {'name': 'tRNA ligase', 
                      'complexes': ['ALAS-MONOMER', 'ARGS-MONOMER', 'ASPS-MONOMER', 
                                    'ASNS-MONOMER', 'CYSS-MONOMER', 'GLURS-MONOMER',
                                    'GLNS-MONOMER', 'GLYS-CPLX', 'HISS-MONOMER', 
                                    'ILES-MONOMER', 'LEUS-MONOMER', 'LYSS-MONOMER',
                                    'METG-MONOMER', 'PHES-CPLX', 'PROS-MONOMER', 
                                    'SERS-MONOMER', 'THRS-MONOMER', 'TRPS-MONOMER',
                                    'TYRS-MONOMER', 'VALS-MONOMER'],
                    'rate_per_sec': 100,
                    'units': 'AA/s',
                    'method':'avg'},
            'glucose_tport': {'name': 'Glucose Transporters', 
                              'complexes': ['CPLX-165', 'CPLX-157'],
                              'rate_per_sec': 200,
                              'units': 'glucose/s',
                              'method': 'sum'},
            'ribosome': {'name': 'Ribosome (50S + 30S)',
                         'complexes': ['CPLX0-3964'],
                         'rate_per_sec': 20,
                         'units': 'AA/s',
                         'method':'sum'},
            'eftu': {'name': 'Elongation Factor EF-Tu', 
                     'complexes': ['EG11037-MONOMER', 'EG11036-MONOMER'],
                     'method': 'sum',
                     'rate_per_sec': 20,
                     'units':'peptide bonds/s'},
            'atp_synthase': {'name': 'F1-F0 ATP Synthase',
                            'complexes':  ['ATPSYN-CPLX'],
                            'method': 'sum', 
                            'rate_per_sec': 100,
                            'units':'atp/s'}, 
            'ndhI': {'name': 'NADH Dehydrogenase I', 
                    'complexes': ['NQOR-CPLX'],
                    'method':'sum',
                    'rate': 4E3,
                    'units': 'protons/s'},
            'fas': {'name': 'Fatty Acid Synthesis',
                    'complexes': ['FABB-MONOMER', 'EG12606', 'EG10277'],
                    'method': 'sum', 
                    'rate': 1,
                    'units':'lipid/s'}}

# %%
complex_df = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr'])):
    for k, v in complexes.items():
        _d = d[d['complex'].isin(v['complexes'])]
        if len(_d) > 0:
            _d = _d.groupby(['complex'])['n_units'].mean().reset_index()
            if v['method'] == 'sum':
                units = _d['n_units'].sum()
                _method = 'sum total'
            if v['method'] == 'avg':
                units = _d['n_units'].mean()
                _method = 'average'

            # assemble a dictionary 
            _data = {'dataset':g[0], 'dataset_name':g[1],
                     'condition':g[2], 'growth_rate_hr':g[3],
                     'n_complex':units, 
                     'rate': v['rate_per_sec'],
                     'rate_units': v['units'],
                     'components':v['complexes'],
                     'shorthand': k,
                     'name': v['name'],
                     'aggregation_method': _method}
        
            complex_df = complex_df.append(_data, ignore_index=True)


complex_df.to_csv('../../data/compiled_estimate_categories.csv', index=False)
# %%
