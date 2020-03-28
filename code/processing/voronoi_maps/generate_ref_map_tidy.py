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
data_group = df_all.groupby(['dataset', 'condition', 'growth_rate_hr'])

#####################################
# Treemap structure and reference map (to provide positional 'seeds' for
# the Voronoi cells).

tree_structure = ['cog_class', 'cog_category', 'gene_name']

cog_dict = dict(zip(df_all.cog_class.unique(),
                    np.arange(len(df_all.cog_class.unique()))))


proteomap_df_ref = \
    gpd.read_file("../../../data/voronoi_map_data/treemap_li_2014_MOPS complete_1.93_ref.geojson",
                                     driver='GeoJSON')

# proteomap_df_ref = \
#     gpd.read_file("../../../data/voronoi_map_data/treemap_schmidt_2016_glucose_0.58.geojson",
#                                      driver='GeoJSON')



count_ = 0
#####################################
# Generate weighted Voronoi map for each dataset
for dets, data in tqdm.tqdm(data_group, desc='Iterating through datasets.'):
    if count_ ==1:
        continue
    count_ += 1

    # only deal with genes > 0.1% of total.
    df_temp = pd.DataFrame()
    for gene, d in data.groupby('gene_name'):
        frac_mass = d.fg_per_cell.sum()/ data.fg_per_cell.sum()
        data_dict = {'frac_mass_tot' : frac_mass,
                        'gene_name' : gene}
        df_temp = df_temp.append(data_dict, ignore_index=True)
    data = data.merge(df_temp, on = 'gene_name')
    data = data[data.frac_mass_tot  >= 0.0001]

    if [dets[0], dets[1]] == ['schmidt_2016', 'glucose']:
        continue

    # initialize DataFrame to store map.
    proteomap_df = pd.DataFrame()

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
                    # S = map.random_points_within(border,  sample_count)
                    S = map.S_find_centroid(proteomap_df_ref[proteomap_df_ref.level == tree_i], tree)

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
            # if tree_i >= 1:
            continue

            area_whole = map.border_map().area
            # print(area_whole)
            for cat, cell in proteomap_df[proteomap_df['level'] == tree_i-1].groupby(tree_structure[tree_i-1]):
                print(cat)
                # if cat != 'RNA processing and modification':
                #     continue
                # For 'Not Assigned' go straigh to gene names
                if cat == 'Not Assigned':
                    tree = 'gene_name'
                    tree_i = 2

                # print('tree level: ', tree_i, ' : ', cat)
                # data to consider for cell
                data_tree = data[data[tree_structure[(tree_i-1)]] == cat]
                print(data_tree)
                border = cell.geometry[cell.index[0]]

                # compute desired cell weightings
                frac = pd.DataFrame()
                for cat_, d in data_tree.groupby(tree):
                    frac_ = d.fg_per_cell.sum()/ data_tree.fg_per_cell.sum()
                    frac_tot = d.fg_per_cell.sum()/ data.fg_per_cell.sum()
                    frac = frac.append(map.data_for_tree(frac_, frac_tot, d, tree_i, cog_dict),
                                ignore_index=True )

                # Set desired cell weighting
                weights = frac.mass_frac.values

                # Set index for each cell
                frac['index'] = np.arange(len(weights))

                # parameters for iteration; try up to 15 times and take best mapping
                error_calc = np.inf
                count = 0

                while (error_calc >= 0.1) and (count <=5):
                    count += 1
                    print(count)
                    try:
                        # number of cells
                        sample_count = len(frac)
                        S = map.random_points_within(border,  sample_count)

                        # Initialize random set of weightings for Voronoi cells.
                        W = (.8 * np.random.random(sample_count) + .2) * (border.area / map.border_map().area) / 1000.0

                        # generate Voronoi cells
                        V = map.compute_power_voronoi_map(S, W, border, 1E-7)

                        # 30 iterations is usually plenty to reach minimum
                        for step in np.arange(30):

                            S, W = map.AdaptPositionsWeights(S, V, W)

                            V = map.compute_power_voronoi_map(S, W, border, 1E-7)

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

                # reset index after handling 'Not Assigned'
                if cat == 'Not Assigned':
                    tree = 'cog_category'
                    tree_i = 1


    # ########################################
    # #  Save to file
    # if proteomap_df.empty:
    #     print('Map not found for: ' + dets[0] + ',' + dets[1] + ',' + str(round(dets[2], 2)))
    # else:
    #     print('count_:', count_, '. Saving map for: ' + dets[0] + ',' + dets[1] + ',' + str(round(dets[2], 2)))
    #     proteomap_df.to_file('../../../data/voronoi_map_data/treemap_' +
    #                 dets[0] + '_' + dets[1] + '_' + str(round(dets[2],2)) + '_ref.geojson',
    #                     driver='GeoJSON')
