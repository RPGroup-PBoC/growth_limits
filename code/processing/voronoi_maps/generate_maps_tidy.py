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

def data_for_tree(frac_, d, tree_i, cog_dict):
    ''' Loads in data details for current mapping (need to make more elegant...)
    '''
    if tree_i == 0:
        data_dict = {'mass_frac' : frac_,
                    'cog_class' : d.cog_class.unique()[0],
                    'cog_category' : None,
                    'cog_dict' : cog_dict[d.cog_class.unique()[0]],
                    'gene_name' : None}
    elif tree_i == 1:
        data_dict = {'mass_frac': frac_,
                    'cog_class' : d.cog_class.unique()[0],
                    'cog_category' : d.cog_category.unique()[0],
                    'cog_dict' : cog_dict[d.cog_class.unique()[0]],
                    'gene_name' : None}
    elif tree_i == 2:
        data_dict = {'mass_frac': frac_,
                   'cog_class' : d.cog_class.unique()[0],
                   'cog_category' : d.cog_category.unique()[0],
                   'cog_dict' : cog_dict[d.cog_class.unique()[0]],
                   'gene_name' : d.gene_name.unique()[0]}

    return data_dict


proteomap_df_ref = \
    gpd.read_file("../../../data/voronoi_map_data/treemap_schmidt_2016_glucose_0.58.geojson",
                                     driver='GeoJSON')


#####################################
# Generate weighted Voronoi map for each dataset
for dets, data in tqdm.tqdm(data_group, desc='Iterating through datasets.'):

    if [dets[0], dets[1]] == ['schmidt_2016', 'glucose']:
        continue

    # initialize DataFrame to store map.
    proteomap_df = pd.DataFrame()

    for tree_i, tree in enumerate(tree_structure):

        if tree_i == 0:
            print('Initialize map.')
            # define boundary of entire map
            border = map.border_map()

            # compute desired cell weightings
            frac = pd.DataFrame()
            for _, d in data.groupby(tree_structure[tree_i]):
                frac_ = d.fg_per_cell.sum()/ data.fg_per_cell.sum()
                frac = frac.append(data_for_tree(frac_, d, tree_i, cog_dict),
                                ignore_index=True )

            # Set index for each cell
            frac['index'] = np.arange(len(frac))

            # Set desired cell weighting
            weights = frac.mass_frac.values

            # parameters for iteration; try up to 15 times and take best mapping
            error_calc = np.inf
            count = 0
            while (error_calc >= 0.1) and (count <=5):
                count += 1
                try:
                    # number of cells
                    sample_count = len(frac)

                    # Set Voronoi points
                    S = map.S_find_centroid(proteomap_df_ref[proteomap_df_ref.level == tree_i], tree)

                    # Initialize random set of weightings for Voronoi cells.
                    W = (.8 * np.random.random(sample_count) + .2)

                    # generate Voronoi cells
                    V, S = map.map_iterator(S, W,  border, weights)

                    # calculate discrepency between desired cell area and
                    # current cell area.
                    A_diff = 0
                    for i, s in enumerate(S):
                        # identify Voronoi cell associated with s
                        for cell_ in V:
                            if (Point(s).within(cell_)):
                                cell =  cell_
                        A_diff += abs(cell.area - border.area * weights[i])

                    if A_diff/(2*border.area) < error_calc:
                        V_cell = V
                        S_final = S

                    # update error value
                    error_calc =  A_diff/(2*border.area)

                except:
                    pass

            # if no errors encountered and map generated, append to DataFrame
            if error_calc != np.inf:
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
            continue
            area_whole = map.border_map().area
            # print(area_whole)
            for cat, cell in proteomap_df[proteomap_df['level'] == tree_i-1].groupby(tree_structure[tree_i-1]):

                # For 'Not Assigned' go straigh to gene names
                if cat == 'Not Assigned':
                    continue
                    # tree = 'gene_name'
                    # tree_i = 2
                print(tree, ' : ',  cat, tree_i)
                # data to consider for cell
                data_tree = data[data[tree_structure[(tree_i-1)]] == cat]

                border = cell.geometry[cell.index[0]]

                # compute desired cell weightings
                frac = pd.DataFrame()
                for cat_, d in data_tree.groupby(tree_structure[(tree_i)]):
                    frac_ = d.fg_per_cell.sum()/ data_tree.fg_per_cell.sum()
                    frac = frac.append(data_for_tree(frac_, d, tree_i, cog_dict),
                                ignore_index=True )

                # Don't attempt to make cells for very low abundant items
                frac = frac[frac.mass_frac >= 0.01]

                # Set index for each cell
                frac['index'] = np.arange(len(frac))

                # Set desired cell weighting
                weights = frac.mass_frac.values

                # parameters for iteration; try up to 15 times and take best mapping
                error_calc = np.inf
                count = 0

                while (error_calc >= 0.1) and (count <=15):
                    count += 1
                    print(count)
                    try:
                        # number of cells
                        sample_count = len(frac)

                        # Set Voronoi points
                        if tree_i == 2:
                            S = map.random_points_within(border,  sample_count)
                        else:
                            map_ref = proteomap_df_ref[proteomap_df_ref['level']==tree_i]
                            map_ref = map_ref[map_ref[tree_structure[tree_i-1]]==cat]
                            S = map.S_find_centroid(map_ref,  tree)

                        # Initialize random set of weightings for Voronoi cells.
                        W = (.8 * np.random.random(sample_count) + .2) * (border.area / area_whole)

                        # generate Voronoi cells
                        V, S = map.map_iterator(S, W,  border, weights)

                        # calculate discrepency between desired cell area and
                        # current cell area.
                        A_diff = 0
                        for i, s in enumerate(S):
                            # identify Voronoi cell associated with s
                            for cell_ in V:
                                if (Point(s).within(cell_)):
                                    cell =  cell_
                            A_diff += abs(cell.area - border.area * weights[i])

                        if A_diff/(2*border.area) < error_calc:
                            V_cell = V
                            S_final = S

                        # update error value
                        error_calc =  A_diff/(2*border.area)


                    except:
                        pass

                # if no errors encountered and map generated, append to DataFrame
                if error_calc != np.inf:
                    gdf = pd.concat([gpd.GeoSeries(cell_) for s in S_final
                                         for cell_ in V_cell
                                         if Point(s).within(cell_)
                                        ]).pipe(gpd.GeoDataFrame)
                    gdf = gdf.rename(columns = {0:'geometry'})
                    gdf['index'] = np.arange(len(frac))
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
        print('saving map for: ' + dets[0] + ',' + dets[1] + ',' + str(round(dets[2], 2)))
        proteomap_df.to_file('../../../data/voronoi_map_data/treemap_' +
                    dets[0] + '_' + dets[1] + '_' + str(round(dets[2],2)) + '.geojson',
                        driver='GeoJSON')
