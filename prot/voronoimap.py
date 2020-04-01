import numpy as np
import pandas as pd
import itertools
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

from shapely.geometry import LineString, MultiLineString, MultiPoint, Point
from shapely.geometry import Polygon, box, MultiPolygon
from shapely.ops import nearest_points, linemerge, unary_union, polygonize
from shapely import affinity

import geopandas as gpd
import fiona
import json


# note that the functions norm2, normalized, get_power_triangulation, and
# get_voronoi_cells were modified from
# https://gist.github.com/marmakoide/45d5389252683ae09c2df49d0548a627#file-laguerre-voronoi-2d-py

def norm2(X):
    """
        Computes

        Parameters
        ----------

        Returns
        -------

    """
    return np.sqrt(np.sum(X ** 2))

def normalized(X):
    return X / norm2(X)

def get_power_triangulation(S, R):
    # Compute the lifted weighted points -NB: THIS IS THE IMPORTANT PART for the weighting
    S_norm = np.sum(S ** 2, axis = 1) - np.array(R)# ** 2

    S_lifted = np.concatenate([S, S_norm[:,None]], axis = 1)

    # Compute the convex hull of the lifted weighted points
    hull = ConvexHull(S_lifted)

    # Extract the Delaunay triangulation from the lower hull
    tri_list = tuple([a, b, c] if is_ccw_triangle(S[a], S[b], S[c]) else [a, c, b]
                    for (a, b, c), eq in zip(hull.simplices, hull.equations)
                    if eq[2] <= 0)
    tri_list_ = tuple([a, b, c] if is_ccw_triangle(S[a], S[b], S[c]) else [a, c, b]
                    for (a, b, c), eq in zip(hull.simplices, hull.equations)
                    if eq[2] >= 0)

    return tri_list, tri_list_


# --- Delaunay triangulation --------------------------------------------------

def get_triangle_normal(A, B, C):
    return normalized(np.cross(A, B) + np.cross(B, C) + np.cross(C, A))


def get_power_circumcenter(A, B, C):
    N = get_triangle_normal(A, B, C)
    return (-.5 / N[2]) * N[:2]

def is_ccw_triangle(A, B, C):
    M = np.concatenate([np.stack([A, B, C]), np.ones((3, 1))], axis = 1)
    return np.linalg.det(M) > 0


# --- Compute Voronoi cells ---------------------------------------------------


def get_voronoi_cells(S, W, tri_list):

    '''
    Compute the segments and half-lines that delimits each Voronoi cell
      * The segments are oriented so that they are in CCW order
      * Each cell is a list of (i, j), (A, U, tmin, tmax) where
         * i, j are the indices of two ends of the segment. Segments end points are
           the circumcenters. If i or j is set to None, then it's an infinite end
         * A is the origin of the segment
         * U is the direction of the segment, as a unit vector
         * tmin is the parameter for the left end of the segment. Can be None, for minus infinity
         * tmax is the parameter for the right end of the segment. Can be None, for infinity
         * Therefore, the endpoints are [A + tmin * U, A + tmax * U]
    '''
    S_norm = np.sum(S ** 2, axis = 1) - np.array(W)
    S_lifted = np.concatenate([S, S_norm[:,None]], axis = 1)
    # Compute the Voronoi points
    V = np.array([get_power_circumcenter(*S_lifted[tri]) for tri in tri_list])
    # Keep track of which circles are included in the triangulation
    vertices_set = frozenset(itertools.chain(*tri_list))

    # Keep track of which edge separate which triangles
    edge_map = { }
    for i, tri in enumerate(tri_list):
        for edge in itertools.combinations(tri, 2):
            edge = tuple(sorted(edge))
            if edge in edge_map:
                edge_map[edge].append(i)
            else:
                edge_map[edge] = [i]

    # For each triangle
    voronoi_cell_map = { i : [] for i in vertices_set }

    for i, (a, b, c) in enumerate(tri_list):
        # For each edge of the triangle
        for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
        # Finite Voronoi edge
            edge = tuple(sorted((u, v)))
            if len(edge_map[edge]) == 2:
                j, k = edge_map[edge]
                if k == i:
                    j, k = k, j

                # Compute the segment parameters
                U = V[k] - V[j]
                U_norm = norm2(U)

                # Add the segment
                voronoi_cell_map[u].append(((j, k), (V[j], U / U_norm, 0, U_norm)))
            else:
            # Infinite Voronoi edge
                # Compute the segment parameters
                A, B, C, D = S[u], S[v], S[w], V[i]
                U = normalized(B - A)
                I = A + np.dot(D - A, U) * U
                W = normalized(I - D)
                if np.dot(W, I - C) < 0:
                    W = -W

                # Add the segment
                voronoi_cell_map[u].append(((edge_map[edge][0], np.inf), (D,  W, 0, np.inf)))
                voronoi_cell_map[v].append(((np.inf, edge_map[edge][0]), (D, -W, np.inf, 0)))



    # Order the segments
    def order_segment_list(segment_list):
        # Pick the first element
        first = min((seg[0][0], i) for i, seg in enumerate(segment_list))[1]

        # In-place ordering
        segment_list[0], segment_list[first] = segment_list[first], segment_list[0]
        for i in range(len(segment_list) - 1):
            for j in range(i + 1, len(segment_list)):
                if segment_list[i][0][1] == segment_list[j][0][0]:
                    segment_list[i+1], segment_list[j] = segment_list[j], segment_list[i+1]
                    break

        return segment_list

    return { i : order_segment_list(segment_list) for i, segment_list in voronoi_cell_map.items() }


def compute_power_voronoi_map(S, W, border, eps):
    '''
    input:
    S - Voronoi points/sites
    W - positive real weights associated with Voronoi points
    eps - for shapely conversion of lines to polygon - needed?
    border - bounds of desired map

    output:
    tri_list: Delaunay triangulation from the lower hull
    V_cell: Set of polygons forming all Voronoi cells

    '''

    # Compute the power triangulation of the circles
    if S.shape[0] > 3:
        tri_list, tri_list_ = get_power_triangulation(S, W)

        # Compute the Voronoi cells
        voronoi_cell_map = get_voronoi_cells(S, W, tri_list)

        # get the lines associated with all the Voronoi cells
        edge_map = { }
        for segment_list in voronoi_cell_map.values():
            for edge, (A, U, tmin, tmax) in segment_list:
                edge = tuple(sorted(edge))
                if edge not in edge_map:
                    if tmax is np.inf:
                        tmax = 20
                    if tmin is np.inf:
                        tmin = -20

                    edge_map[edge] = (A + tmin * U, A + tmax * U)

        V_map = MultiLineString(list(edge_map.values())).buffer(eps)

        V_cell = border.difference(V_map)

        while len(S) != len(V_cell):
            V_cell = VoronoiMap_fix(S, V_cell)


    elif S.shape[0] <= 3:
        S_points = MultiPoint(S[:2])
        bisect = affinity.rotate(LineString(S_points), 90, 'center')
        bisect = affinity.scale(bisect, xfact=100, yfact=100)

        borders = unary_union(linemerge([border.boundary, bisect]))
        V_cell = MultiPolygon(polygonize(borders))

        while len(S) != len(V_cell):
            V_cell = VoronoiMap_fix(S, V_cell)

    return V_cell


def VoronoiMap_fix(S, V_cell):
    '''
    Use to shift Voronoi sites and adapt weights W
    input:
    V_cell - Set of polygons forming all Voronoi cells
    S: Voronoi points

    output:
    V_cell_corr - Corrected Voronoi polygons
    '''

    S_points = MultiPoint(S)

    V_cell_ = MultiPolygon([cell for cell in V_cell if ((cell.intersection(S_points)).type  != 'MultiPoint')])
    V_ = MultiPolygon([cell for cell in V_cell if ((cell.intersection(S_points)).type  == 'MultiPoint')])

    while len(V_) >=1:
        for cell in V_:
            V_ = MultiPolygon([P for P in V_ if P != cell])
            # identify the multipe voronoi points in cell
            P = cell.intersection(S_points)
            P_coords = [((p_.xy[0][0], p_.xy[1][0])) for p_ in P]
            for s in P:
                # find point closest to s and split cell using line orthogonal to
                # their center line
                s_coords = ((s.xy[0][0], s.xy[1][0]))
                s_other = [tuple(s_) for s_ in P_coords if tuple(s_) != tuple(s_coords)]
                NN = nearest_points((s),MultiPoint(s_other))
                bisect = affinity.rotate(LineString(NN), 90, 'center')
                bisect = affinity.scale(bisect, xfact=100, yfact=100)

                borders = unary_union(linemerge([cell.boundary, bisect]))
                cell_split = polygonize(borders)

                for cell_ in cell_split:
                    if ((cell_.intersection(S_points)).type  == 'MultiPoint'):
                        poly=[]
                        for pol in V_:
                            poly.append(pol)
                        poly.append(cell_)
                        V_ = MultiPolygon(poly)
                    else:
                        poly=[]
                        for pol in V_cell_:
                            poly.append(pol)
                        poly.append(cell_)
                        V_cell_ = MultiPolygon(poly)
                break
    return V_cell_


def AdaptPositionsWeights(S, V_cell, W):
    '''
    Use to shift Voronoi sites and adapt weights W
    input:
    V_cell - Set of polygons forming all Voronoi cells
    W - positive real weights associated with Voronoi points

    output:
    S_new: updated Voronoi points
    W_new: updates weights
    '''

    S_ = []
    dist_border = []
    for s in S:
        # identify Voronoi cell associated with s
        for cell_ in V_cell:
            if (Point(s).within(cell_)):
                cell =  cell_
        # Move s to centroid of cell
        for pp in cell.centroid.coords:
            S_.append(pp)
#             break
            dist_border.append(abs(cell.exterior.distance(Point(pp))))
    # S_new = np.array(S_)
    S_new = np.array(list(map(list, S_)))
    W_new = [np.min([np.sqrt(W[i]),dist_border[i]])**2 for i in np.arange(len(W))]

    return S_new, W_new



def AdaptWeights(V_cell, S, bound, W, w_desired, err = 0.001):
    '''
    Use to adapt weights W
    input:
    V_cell - Set of polygons forming all Voronoi cells
    S - Voronoi points
    W - positive real weights associated with Voronoi points
    w_desired - desired weighting of each cell
    err - error threshold for values in W; need to understand better.

    output:
    W_new: updates weights
    '''
    W_new = []

    for i, s in enumerate(S):
        s_other = [tuple(s_) for s_ in S if tuple(s_) != tuple(s)]
        NN = nearest_points(Point(s),MultiPoint(s_other))

        for cell_ in V_cell:
            if (Point(s).within(cell_)):
                cell = cell_

        A_current = cell.area
        A_desired = bound.area * w_desired[i]

        f_adapt = A_desired / A_current

        w_new_ = np.sqrt(W[i]) * f_adapt

        w_max = abs(LineString(NN).length)

        W_ = np.min([w_new_, w_max])**2

        W_new.append(np.max([W_, err]))

    return W_new

# def compute_power_voronoi_map_(S, W, border, eps = 1E-6):
#     result = None
#     while result is None:
#         try:
#             V_cell = compute_power_voronoi_map(S, W, border, eps = 1E-6)
#             result = True
#             return V_cell
#
#         except:
#             pass

def map_iterator(S, W,  border, weights):
    # generate initial Voronoi cells
    # print(S)
    V = compute_power_voronoi_map(S, W, border, 1E-7)

    # 30 iterations is usually plenty to reach minimum
    for i in np.arange(30):
        S, W = AdaptPositionsWeights(S, V, W)

        V = compute_power_voronoi_map(S, W, border, 1E-7)

        W = AdaptWeights(V, S, border, W, weights, 1E-7)

        V = compute_power_voronoi_map(S, W, border, 1E-7)
        # print(V)
    return V, S

def border_map():
    b_s = 15
    return box(-b_s, -b_s, b_s, b_s)

# def map_iterator(weights, border):
#     # number of cells
#     sample_count = len(weights)
#
#     # initialize random Voronoi sites and weightings
#     # S contains circles center, R = sqrt(W) contains
#     # circles radius for weighting the Voronoi cells
#     S_ = random_points_within(border, sample_count)
#     S = np.array([([s.x, s.y]) for s in S_])
#
#     W = .8 * np.random.random(sample_count) + .2
#
#     # generate initial Voronoi cells
#     V_cell = compute_power_voronoi_map_(S, W, border, eps = 1E-6)
#
#     # Perform iteration on initialized V_cell
#     error_diff = np.inf
#     error_calc = 100.0
#     count = 0
# #     while 0 > error_diff >= 0.00001:
#     while np.logical_and(error_calc >= 0.08, count <=101):
#         print(np.logical_and(error_calc >= 0.08, count <=101))
#         print(error_diff, error_calc, count)
#         S, W = AdaptPositionsWeights(S, V_cell, W)
#
#         V_cell = compute_power_voronoi_map_(S, W, border, eps = 1E-6)
#
#         W = AdaptWeights(V_cell, S, border, W, weights, err = 1E-5)
#
#         V_cell = compute_power_voronoi_map_(S, W, border, eps = 1E-6)
#
#         # calculate discrepency between desired cell area and
#         # current cell area.
#         A_diff = 0
#         for i,  s in enumerate(S):
#             # identify Voronoi cell associated with s
#             for cell_ in V_cell:
#                 if (Point(s).within(cell_)):
#                     cell =  cell_
#             A_diff += abs(cell.area - border.area * weights[i])
#
#         error_diff = error_calc -  A_diff/(2*border.area)
#         error_calc =  A_diff/(2*border.area)
#         count += 1
#
#
#
#     return V_cell

# selection of random initialization points
def random_points_within(poly, num_points):
    min_x, min_y, max_x, max_y = poly.bounds

    points = []

    while len(points) < num_points:
        x, y = [np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)]
        random_point = Point([x, y])
        if (random_point.within(poly)):
            points.append([x, y])

    S = np.array(list(map(list, points)))
    return S


def S_find_centroid(proteomap, tree, tree_list, border, random_shift = False):
    """
    Generates the Voronoi positions based on  the centroids of cells in
    a  reference map.

    Parameters
    ----------
    proteomap: pandas DataFrame
        Dataframe with Voronoi map information.
    treeID: str
        Key by which to indicate the tree level category (e.g. 'cog_class').
    level: int
        Key to indicate which level of the map to consider.

    Returns
    -------
    S: (N,  2) numpy.ndarray
        Array of N Voronoi points.
    """

    S_ = []

    for cell_id, d in proteomap.groupby(tree):
        if cell_id in tree_list:
            S_.append([d.geometry.centroid.x.values[0], d.geometry.centroid.y.values[0]])

    num_missing =  len(tree_list) - len(S_)
    if num_missing >= 1:
        S_.append(random_points_within(border,  num_missing))

    if random_shift == True:
        mu, sigma = 1.0, 0.35 # mean and standard deviation
        S_ = [s*np.random.normal(mu, sigma, 2) for s in S_]

    S = np.array(list(map(list, S_)))

    return S


def S_transform(S_shapely, border, translation):
    """
    Performs translation and scaling of Voronoi positions based on the difference
    of current cell and reference cell.

    Parameters
    ----------
    S_shapely: MultiPoint shapely type

    Returns
    -------
    S: (N,  2) numpy.ndarray
        Array of N Voronoi points.
    """
    if S_shapely.within(border) == False:
        S_shapely = shapely.affinity.translate(S_shapely, xoff=translation[0], yoff=translation[1])

        while S_shapely.within(border) == False:
            S_shapely = shapely.affinity.scale(S_shapely, xfact=0.9, yfact=0.9)

        S = np.array([[s.x,s.y] for s in S_shapely])

        return S
    else:
        S = np.array([[s.x,s.y] for s in S_shapely])
        return S


def data_for_tree(frac, frac_tot, d, tree_i, cog_dict):
    ''' Loads in data details for current mapping (need to make more elegant...)
    '''
    if tree_i == 0:
        data_dict = {'mass_frac' : frac,
                    'mass_frac_tot' : frac_tot,
                    'cog_class' : d.cog_class.unique()[0],
                    'cog_category' : None,
                    'cog_dict' : cog_dict[d.cog_class.unique()[0]],
                    'gene_name' : None}
    elif tree_i == 1:
        data_dict = {'mass_frac': frac,
                    'mass_frac_tot' : frac_tot,
                    'cog_class' : d.cog_class.unique()[0],
                    'cog_category' : d.cog_category.unique()[0],
                    'cog_dict' : cog_dict[d.cog_class.unique()[0]],
                    'gene_name' : None}
    elif tree_i == 2:
        data_dict = {'mass_frac': frac,
                    'mass_frac_tot' : frac_tot,
                   'cog_class' : d.cog_class.unique()[0],
                   'cog_category' : d.cog_category.unique()[0],
                   'cog_dict' : cog_dict[d.cog_class.unique()[0]],
                   'gene_name' : d.gene_name.unique()[0]}

    return data_dict
