"""
This script takes a dataframe of protein complexes obtained using Ecocyc's
SmartTable tool, converts it into a .csv file, and uses the Object IDs to
scrape EcoCyc to identify the subunit numbers.

TODO: ensure it is robust; double check numbers
"""
#%%
import numpy as np
import pandas as pd
import tqdm

# tools for website scraping
import requests
from bs4 import BeautifulSoup
import re

# Grab all enzyme/ protein complexes in Ecocyc.
df_complex = pd.read_csv('../../../data/Ecocyc_20200205_All_protein_complexes_of_E._coli_K-12_substr._MG1655.txt', delimiter='	')

# clean up and save as .csv
df_complex_ = pd.DataFrame()

for item, data in tqdm.tqdm(df_complex.groupby('Genes'), desc='Cleaning up Ecocyc data'):
##     if len(data.Genes.values[0].split('//'))>1:
    for gene in  data.Genes.values[0].split('//'):
        complex_name = data.Complexes_Proteins.values[0]

        # clean up the complex name from html to latex
        if '<i>' in complex_name:
            complex_name = complex_name.replace('<i>', '{\it ')
            complex_name = complex_name.replace('</i>', '}')
        if '<sup>' in complex_name:
            complex_name = complex_name.replace('<sup>', '$^')
            complex_name = complex_name.replace('</sup>', '$')
        if '<sub>' in complex_name:
            complex_name = complex_name.replace('<sub>', '$_{')
            complex_name = complex_name.replace('</sub>', '}$')
        if '&' in complex_name:
            complex_name = complex_name.replace('&', '$\\')
            complex_name = complex_name.replace(';', '$')
        if '<small>' in complex_name:
            complex_name = complex_name.replace('<small>', '{\small ')
            complex_name = complex_name.replace('</small>', '}')

        data_list = {'complex' : complex_name,
               'gene' : re.sub('[\W_]+', '', gene),
               'object_id': data.Object_ID.values[0]}
        df_complex_ = df_complex_.append(data_list,  ignore_index=True)

df_complex_.to_csv('../../../data/Ecocyc_20200205_All_protein_complexes_of_E._coli_K-12_substr_MG1655.csv')


# Website scraping
df_subunits = pd.DataFrame()
for object_id in tqdm.tqdm(df_complex_.object_id.unique(), desc='Iterating through complexes'):
    URL = 'https://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=' + object_id
    page = requests.get(URL)

    soup = BeautifulSoup(page.content, 'html.parser')

    # mark status (0 = subunit details not found; 1 = found)
    check = 0

    # attempt to find details in 'POLYPEPTIDE' class
    for string in soup.find_all(class_='POLYPEPTIDE'):
        if object_id in repr(string):
            ## get formula
            full_seq = repr(string).split('=  ')[1].split(']')[:-1]
            full_seq = re.split("(])",repr(string).split('=  ')[1])
            full_seq = ''.join(full_seq)
            full_seq = ''.join(full_seq.split('\n')[:-1])
            ####
            for obj in repr(string).split('=  ')[1].split('[')[1:]:

                if re.findall(r'\d+', obj) == []:
                    subunit_count = 1
                else:
                    subunit_count = int(re.findall(r'\d+', obj)[0])

                data_list = {'object_id' : object_id,
                          'gene' : obj.split(']')[0][0].lower() + obj.split(']')[0][1:],
                          'subunit_count' : subunit_count,
                          'formula' : full_seq}
                df_subunits = df_subunits.append(data_list,ignore_index=True)
                check = 1

    if check == 0:
        # attempt to find details in 'ENZYME' class
        for string in soup.find_all(class_='ENZYME'):
            if object_id in repr(string):
                if '=  ' in repr(string):
                    ## get formula
                    full_seq = repr(string).split('=  ')[1].split(']')[:-1]
                    full_seq = re.split("(])",repr(string).split('=  ')[1])
                    full_seq = ''.join(full_seq)
                    full_seq = ''.join(full_seq.split('\n')[:-1])
                    ####
                    for obj in repr(string).split('=  ')[1].split('[')[1:]:
#                             print(obj)
                        if re.findall(r'\d+', obj) == []:
                            subunit_count = 1
                        else:
                            subunit_count = int(re.findall(r'\d+', obj)[0])
                        data_list = {'object_id' : object_id,
                                  'gene' : obj.split(']')[0][0].lower() + obj.split(']')[0][1:],
                                  'subunit_count' : subunit_count ,
                                  'formula' : full_seq}
                        df_subunits = df_subunits.append(data_list,ignore_index=True)
                        check = 1
    if check == 0:
        # attempt to find details in ECOLI general orgid from object_ID
        URL = 'https://ecocyc.org/gene?orgid=ECOLI&id=' + object_id
        for gene in df_complex_[df_complex_.object_id == object_id].gene.unique():
            check = 0
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            test = soup.findAll('a')
            for items in test:

                if '['+gene[0].upper()+gene[1:]+']' in repr(items):
                    content = repr(items).split(';/b&gt;      ')
                    for sec in content:
                        if '['+gene[0].upper()+gene[1:]+']' in sec:
                            ## get formula
                            full_seq = '[' + ''.join(sec.split('[')[1:])
                            full_seq = ''.join(full_seq.split(', WIDTH')[0][:-1])
                            ####
                            test = repr(items).split('['+gene[0].upper()+gene[1:]+']')[1][:12]

                            if re.findall(r'\d+', test.split(',')[0]) == []:
                                subunit_count = 1
                            else:############ need to double check
                                subunit_count = int(re.findall(r'\d+', test)[0])
                            data_in = {'object_id' : object_id,
                              'gene' : gene,
                              'subunit_count' : subunit_count,
                              'formula' : full_seq}
                            df_subunits = df_subunits.append(data_in, ignore_index=True)
                            check = 1
                            break

    if check == 0:
        # details not found; set values to np.nan
        print('Did not find object: ', object_id)
        for gene in df_complex_[df_complex_.object_id == object_id].gene.unique():
            data_list = {'object_id' : object_id,
              'gene' : gene,
              'subunit_count' : np.nan}
            df_subunits = df_subunits.append(data_list, ignore_index=True)


# save to file
df_subunits.to_csv('../../../data/protein_complexes_MG1655_subunit_numbers.csv')
