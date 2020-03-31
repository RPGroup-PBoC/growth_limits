import glob
import tqdm

import numpy as np
import pandas as pd

# import itertools
# import scipy
# from scipy.spatial import ConvexHull
# from scipy.spatial import Delaunay

import shapely
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point
from shapely.geometry import Polygon, box, MultiPolygon
from shapely.ops import nearest_points, linemerge, unary_union, polygonize
from shapely import affinity

import geopandas as gpd
import fiona
import json

import prot.voronoimap as map

#####################################
# Load the dataset with aboslute measurements

df_all = pd.read_csv('../../../data/compiled_absolute_measurements.csv')

# use dataset from Li et al. complete media, since this is most complete dataset (right?)
df = df_all[df_all.dataset == 'li_2014']
df = df[df.condition == 'MOPS complete']

# generate a few more subcategories in the 'information storage and processing'
# COG class :
## ------------------------
# transcription
## ------------------------
## 'transcription factors and nucleoid  associated  proteins' : GO:0003677) - DNA binding
## 'RNA polymerase and sigma factors' : GO:0016987; sigma factor activity or
        ## GO:0003899; DNA-directed 5'-3' RNA polymerase activity, DNA-directed RNA polymerase,
## 'transcription related' : else

## ------------------------
# 'translation, ribosomal structure and biogenesis' :
## ------------------------
## 'Ribosome' : GO:0005840; ribosome, Ribonucleoprotein; Ribosomal protein
## 'translation, ribosomal biogenesis' : else

## ------------------------
# Replication, recombination, and repair
## ------------------------
## 'DNA polymerase' : GO:0003887; DNA-directed DNA polymerase activity
## 'DNA replication related' : GO:0006260; DNA replication
## 'DNA Recombination and repair' : else

df_ref = pd.DataFrame()

for cog_cat, _d in df.groupby('cog_category'):
    if cog_cat == 'transcription':
        d = _d[_d['go_terms'].astype(str).str.contains('GO:0016987|GO:0003899')] #GO:0016987; sigma factor activity
        d = d.replace({'transcription': 'RNA polymerase and sigma factors'})
        df_ref = df_ref.append(d, ignore_index = True)

        d = _d[~_d['go_terms'].astype(str).str.contains('GO:0016987|GO:0003899')]
        d = d[d['go_terms'].astype(str).str.contains('GO:0003677')]
        d = d.replace({'transcription': \
                'transcription factors and nucleoid associated proteins'})
        df_ref = df_ref.append(d, ignore_index = True)

        d = _d[~_d['go_terms'].astype(str).str.contains('GO:0016987|GO:0003899')]
        d = d[~d['go_terms'].astype(str).str.contains('GO:0003677')]
        d = d.replace({'transcription': 'transcription related'})
        df_ref = df_ref.append(d, ignore_index = True)

    elif cog_cat == 'translation, ribosomal structure and biogenesis':
        d = _d[_d['go_terms'].astype(str).str.contains('GO:0005840')]
        d = d.replace({'translation, ribosomal structure and biogenesis' : \
                                    'ribosome'})
        df_ref = df_ref.append(d, ignore_index = True)

        d = _d[~_d['go_terms'].astype(str).str.contains('GO:0005840')]
        d = d.replace({'translation, ribosomal structure and biogenesis' : \
                 'translation, ribosomal biogenesis'})
        df_ref = df_ref.append(d, ignore_index = True)

    elif cog_cat == 'Replication, recombination, and repair':
        d = _d[_d['go_terms'].astype(str).str.contains('GO:0003887')]
        d = d.replace({'Replication, recombination, and repair': \
                    'DNA polymerase'})
        df_ref = df_ref.append(d, ignore_index = True)

        d = _d[~_d['go_terms'].astype(str).str.contains('GO:0003887')]
        d = d[d['go_terms'].astype(str).str.contains('GO:0006260')]
        d = d.replace({'Replication, recombination, and repair': \
                    'DNA replication related'})
        df_ref = df_ref.append(d, ignore_index = True)

        d = _d[~_d['go_terms'].astype(str).str.contains('GO:0003887')]
        d = d[~d['go_terms'].astype(str).str.contains('GO:0006260')]
        d = d.replace({'Replication, recombination, and repair': \
                    'DNA Recombination and repair'})
        df_ref = df_ref.append(d, ignore_index = True)

    else:
        df_ref = df_ref.append(_d, ignore_index = True)

print(df_ref.cog_category.unique())


data_group = df_ref.groupby(['dataset', 'condition', 'growth_rate_hr'])

#####################################
# Treemap structure and reference map (to provide positional 'seeds' for
# the Voronoi cells).

tree_structure = ['cog_class', 'cog_category', 'gene_name']

cog_dict = dict(zip(df_all.cog_class.unique(),
                    np.arange(len(df_all.cog_class.unique()))))


proteomap_df = \
    gpd.read_file("../../../data/voronoi_map_data/treemap_li_2014_MOPS complete_1.93_ref.geojson",
                                     driver='GeoJSON')
proteomap_df = proteomap_df[proteomap_df.level <= 1]
# proteomap_df_ref = \
#     gpd.read_file("../../../data/voronoi_map_data/treemap_schmidt_2016_glucose_0.58.geojson",
#                                      driver='GeoJSON')




#####################################
# Generate weighted Voronoi map for each dataset
for dets, data in tqdm.tqdm(data_group, desc='Iterating through datasets.'):

    # Use this frac_mass_tot only for tree_i = 2 (gene_names)
    # only deal with genes > 0.1% of total.
    df_temp = pd.DataFrame()
    for gene, d_ in data.groupby('gene_name'):
        frac_mass = d_.fg_per_cell.sum()/ data.fg_per_cell.sum()
        data_dict = {'frac_mass_tot' : frac_mass,
                        'gene_name' : gene}
        df_temp = df_temp.append(data_dict, ignore_index=True)
    data_genes = data.copy().merge(df_temp, on = 'gene_name')
    # data_genes = data_genes[data_genes.frac_mass_tot  >= 0.0001]

    # # initialize DataFrame to store map.
    # proteomap_df = pd.DataFrame()

    for tree_i, tree in enumerate(tree_structure):

        if tree_i == 0:
            continue
            print('Initialize map.')

            # define boundary of entire map
            border = map.border_map()

            # compute desired cell weightings
            frac = pd.DataFrame()
            for _, d in data.groupby(tree_structure[tree_i]):
                frac_ = d.fg_per_cell.sum()/ data.fg_per_cell.sum()
                frac = frac.append(map.data_for_tree(frac_, frac_, d, tree_i, cog_dict),
                                ignore_index=True )


            # Set index for each cell
            frac['index'] = np.arange(len(frac))

            # Set desired cell weighting
            weights = frac.mass_frac.values

            # parameters for iteration; try up to 15 times and take best mapping
            error_calc = np.inf
            count = 0
            while (error_calc >= 0.1) and (count <=15):
                count += 1
                try:
                    # number of cells
                    sample_count = len(frac)

                    # Set Voronoi points
                    S = map.random_points_within(border,  sample_count)
                    # S = map.S_find_centroid(proteomap_df_ref[proteomap_df_ref.level == tree_i], tree)

                    # Initialize random set of weightings for Voronoi cells.
                    W = (.8 * np.random.random(sample_count) + .2)

                    # generate Voronoi cells
                    V = map.compute_power_voronoi_map(S, W, border, 1E-7)
                    # print(V,S)
                    # 30 iterations is usually plenty to reach minimum
                    for step in np.arange(30):

                        S, W = map.AdaptPositionsWeights(S, V, W)

                        V = map.compute_power_voronoi_map(S, W, border, 1E-7)

                        W = map.AdaptWeights(V, S, border, W, weights, 1E-3)

                        V = map.compute_power_voronoi_map(S, W, border, 1E-7)

                        # calculate discrepency between desired cell area and
                        # current cell area.
                        A_diff = 0
                        for i, s in enumerate(S):
                            # identify Voronoi cell associated with s
                            for cell_ in V:
                                if (Point(s).within(cell_)):
                                    cell =  cell_
                            A_diff += abs(cell.area - border.area * weights[i])
                        print(A_diff/(2*border.area))

                        if A_diff/(2*border.area) < error_calc:
                            V_cell = V
                            S_final = S

                            # update error value
                            error_calc =  A_diff/(2*border.area)

                except:
                    pass

            # if no errors encountered and map generated, append to DataFrame
            if error_calc != np.inf:
                # print(error_calc)
                gdf = pd.concat([gpd.GeoSeries(cell_) for s in S_final
                                     for cell_ in V_cell
                                     if Point(s).within(cell_)
                                    ]).pipe(gpd.GeoDataFrame)
                gdf = gdf.rename(columns = {0:'geometry'})
                gdf['index'] = np.arange(len(frac))
                gdf['level'] = tree_i

                proteomap_df = proteomap_df.append(gdf.merge(frac, on=['index']))
            else:
                print('Did not find map for top level.')
                break

        else:


            area_whole = map.border_map().area
            # print(area_whole)
            for cat, cell in proteomap_df[proteomap_df['level'] == tree_i-1].groupby(tree_structure[tree_i-1]):
                print(cat)
                # if cat != 'RNA processing and modification':
                #     continue
                # For 'Not Assigned' go straigh to gene names
                if tree_i == 1:
                    continue
                    # if cat == 'Not Assigned':
                    #     tree = 'gene_name'
                    #     tree_i = 2
                    # else:
                    #     continue

                # if tree_i == 2:
                #     data = data_genes

                # data to consider for cell
                data_tree = data[data[tree_structure[(tree_i-1)]] == cat]

                border = cell.geometry[cell.index[0]]

                # compute desired cell weightings
                frac = pd.DataFrame()
                for cat_, d in data_tree.groupby(tree):
                    frac_ = d.fg_per_cell.sum()/ data_tree.fg_per_cell.sum()
                    frac_tot = d.fg_per_cell.sum()/ data.fg_per_cell.sum()
                    frac = frac.append(map.data_for_tree(frac_, frac_tot, d, tree_i, cog_dict),
                                ignore_index=True )

                # sort and take top 95% by weight
                frac = frac.sort_values('mass_frac', ascending=False)
                frac['cumsum'] = np.cumsum(frac.mass_frac.values)
                # frac = frac[frac['cumsum'] <= 0.95]

                if frac.empty:
                    continue

                # Set index for each cell
                frac['index'] = np.arange(len(frac))

                # truncate genes due to computational challenge (if needed)
                if frac['index'].max() >= 51:
                    frac = frac[frac['index'] <= 50]

                # Set desired cell weighting
                weights = frac.mass_frac.values
                print(weights)


                # parameters for iteration; try up to 15 times and take best mapping
                error_calc = np.inf
                count = 0

                while (error_calc >= 0.15) and (count <=15):
                    count += 1
                    print(count)
                    try:
                        # number of cells
                        sample_count = len(frac)
                        S = map.random_points_within(border,  sample_count)

                        # Initialize random set of weightings for Voronoi cells.
                        W = (.8 * np.random.random(sample_count) + .2) * (border.area / map.border_map().area)#/ 1000.0

                        # generate Voronoi cells
                        V = map.compute_power_voronoi_map(S, W, border, 1E-7)

                        # 30 iterations is usually plenty to reach minimum
                        for step in np.arange(30):

                            S, W = map.AdaptPositionsWeights(S, V, W)

                            V = map.compute_power_voronoi_map(S, W, border, 1E-7)
                            if cat == 'information storage and processing':
                                W = map.AdaptWeights(V, S, border, W, weights, 1E-5)
                            else:
                                W = map.AdaptWeights(V, S, border, W, weights, 1E-6)

                            V = map.compute_power_voronoi_map(S, W, border, 1E-7)

                            # calculate discrepency between desired cell area and
                            # current cell area.
                            A_diff = 0
                            for i, s in enumerate(S):
                                # identify Voronoi cell associated with s
                                for cell_ in V:
                                    if (Point(s).within(cell_)):
                                        cell =  cell_
                                A_diff += abs(cell.area - border.area * weights[i])
                            print(A_diff/(2*border.area))

                            if A_diff/(2*border.area) < error_calc:
                                V_cell = V
                                S_final = S

                                # update error value
                                error_calc =  A_diff/(2*border.area)


                    except:
                        pass

                # if no errors encountered and map generated, append to DataFrame
                if error_calc != np.inf:
                    print(error_calc)
                    gdf = pd.concat([gpd.GeoSeries(cell_) for s in S_final
                                         for cell_ in V_cell
                                         if Point(s).within(cell_)
                                        ]).pipe(gpd.GeoDataFrame)
                    gdf = gdf.rename(columns = {0:'geometry'})
                    gdf['index'] = np.arange(len(weights))
                    gdf['level'] = tree_i

                    proteomap_df = proteomap_df.append(gdf.merge(frac, on=['index']))
                else:
                    print('Did not find map for:' + cat)

                # # reset index after handling 'Not Assigned'
                # if cat == 'Not Assigned':
                #     tree = 'cog_category'
                #     tree_i = 1


    ########################################
    #  Save to file
    if proteomap_df.empty:
        print('Map not found for: ' + dets[0] + ',' + dets[1] + ',' + str(round(dets[2], 2)))
    else:
        print('. Saving map for: ' + dets[0] + ',' + dets[1] + ',' + str(round(dets[2], 2)))
        proteomap_df.to_file('../../../data/voronoi_map_data/treemap_' +
                    dets[0] + '_' + dets[1] + '_' + str(round(dets[2],2)) + '_ref_.geojson',
                        driver='GeoJSON')
