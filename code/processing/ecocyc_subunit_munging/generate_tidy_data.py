"""
This script takes a dataframe of protein complexes obtained using Ecocyc's
SmartTable tool, converts it into a .csv file, and uses the Object IDs to
scrape EcoCyc to identify the subunit numbers.

TODO: Currently broken (3/11/20 - GC)
"""
#%%
import numpy as np
import pandas as pd
import tqdm
import skimage.measure

# tools for website scraping
import requests
from bs4 import BeautifulSoup
import re

# Grab all enzyme/ protein complexes in Ecocyc.
########################## Website scraping
##################################################
df_complex = pd.read_csv('../../../data/ecocyc_raw_data/Ecocyc_20200205_All_protein_complexes_of_E._coli_K-12_substr._MG1655.txt', delimiter='	')
manual_annotation = pd.read_csv('../../../data/ecocyc_raw_data/manually_annotated_complexes.csv') 
annotations = pd.read_csv('../../../data/ecoli_genelist_master.csv')

#%%
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

#%%
########################## Website scraping
##################################################

##########
# Identify formula associated with each object_id
##########

df_subunits = pd.DataFrame()
missing_cplx = []
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
            if '<sub>' in full_seq:
                full_seq = full_seq.replace('<sub>', '$_{')
                full_seq = full_seq.replace('</sub>', '}$')

            data_list = {'object_id' : object_id,
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
                    if '<sub>' in full_seq:
                        full_seq = full_seq.replace('<sub>', '$_{')
                        full_seq = full_seq.replace('</sub>', '}$')

                    data_list = {'object_id' : object_id,
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
                    content = repr(items).split('&lt;br&gt;')
                    for sec in content:
                        if '['+gene[0].upper()+gene[1:]+']' in sec:

                            ## get formula

                            full_seq = '[' + ''.join(sec.split('[')[1:])
#                             print(full_seq)
                            full_seq = ''.join(full_seq.split(', WIDTH')[0])

                            if '&lt;SUB&gt;' in full_seq:
                                full_seq = full_seq.replace('&lt;SUB&gt;', '$_{')
                                full_seq = full_seq.replace('&lt;/SUB&gt;', '}$')

                            if full_seq[-1] == '\'':
                                full_seq = full_seq[:-1]
                            data_in = {'object_id' : object_id,
                              'formula' : full_seq}
                            df_subunits = df_subunits.append(data_in, ignore_index=True)
                            check = 1
                        if check == 1:
                            break

    if check == 0:
        # details not found; set values to np.nan
        print('Did not find object: ', object_id)
        missing_cplx.append(object_id)
        # for gene in df_complex_[df_complex_.object_id == object_id].gene.unique():
        #     data_list = {'object_id' : object_id,'formula' : ''}
        #     df_subunits = df_subunits.append(data_list, ignore_index=True)


# %%
# Load the manually annotated missing complexes and append.
df_subunits_complete = pd.concat([df_subunits, manual_annotation], sort=False)
df_subunits_complete

#%%
# Use formula to determine subunit counts
df_subunits_nums = pd.DataFrame()
_dfs = []
for obj, data  in tqdm.tqdm(df_subunits_complete.groupby(['object_id', 'formula']), desc='Calculating subunit numbers'):
    obj_genes = df_complex_[df_complex_.object_id==obj[0]].gene.unique()
    gene_chars = ''.join(list(obj_genes)).lower()
    formula_ = data.formula.unique()[0]
    formula_ = formula_.replace('$_{', '')
    formula_ = formula_.replace('}$', '')
    formula_ = formula_.replace(' ', '')
    formula_ = formula_.replace('(', '[')
    formula_ = formula_.replace(')', ']')
    brackets = np.where(np.array([(f == '[') | (f == ']') for f in formula_]) ==1)[0]
    numbers = skimage.measure.label(np.array([f.isnumeric() for f in formula_.lower()]))
    parsed_genes = [formula_[brackets[i] + 1:brackets[i+1]] for i in range(len(brackets) - 1) if len(formula_[brackets[i] +1:brackets[i + 1]])> 0]
    subunits = []
    if (len(numbers) - 1) == 0:
        subunits = [1 for _ in range(len(parsed_genes))]
    _df = pd.DataFrame([]) 
    _df['genes'] = parsed_genes
    _df['n'] = subunits
    _df['formula'] = obj[1]
    _dfs.append(_df)
_df = pd.concat(_dfs, sort=False)
#%%

#
#    for gene in re.split("[^a-zA-Z]*",formula_)[1:-1]:
#        # check that gene identified is actually expected; otherwise skip
#        if gene.lower() not in [name.lower() for name in obj_genes]:
#            continue
#        formula = formula_
#        # remove all other genes in the formula and calculate subunit total count.
#        for item in re.split("[^a-zA-Z]*",formula)[1:-1]:
#            print(item)
#            if item != gene:
#                ind = formula.find(item)
#                if formula[np.max([0,(ind-1)]):(ind+len(item)+1)][-1] == ')':
#                    step = 1
#                else:
#                    step = 0
#                if formula[(ind-1):(ind+len(item)+step+1)] == formula[(ind-1):]:
#                    formula = formula.replace(formula[(ind-1-step):], '')
#
#                elif formula[np.max([0,(ind-1)]):(ind+len(item)+step+2)][-1].isnumeric():
#                    if formula[np.max([0,(ind-1)]):(ind+len(item)+step+3)][-1].isnumeric():
#                        formula = formula.replace(formula[np.max([0,(ind-1-step)]):(ind+len(item)+step+3)], '')
#                    else:
#                        formula = formula.replace(formula[np.max([0,(ind-1-step)]):(ind+len(item)+step+2)], '')
#                else:
#                    formula = formula.replace(formula[np.max([0,(ind-1-step)]):(ind+len(item)+step+2)][:-1], '')
#
#        data_list = {'object_id':  obj[0], 'formula' : obj[1],
#                    'gene' : gene, 'subunit_count': np.prod([ int(x) for x in re.findall(r'[0-9]+', formula)] ) }
#        df_subunits_nums = df_subunits_nums.append(data_list, ignore_index=True)
#
#
##%%
## save to file
#df_subunits_nums.to_csv('../../../data/protein_complexes_MG1655_subunit_numbers.csv')
#