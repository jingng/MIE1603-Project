import copy

import numpy as np
import pandas as pd
import geopandas as gpd
from osmnx import utils_geo as ug
from sklearn.neighbors import BallTree
from shapely import MultiPoint, transform, extract_unique_points, plotting

from runtrack.graph import distance, utils


class Candidate:
    """ Class to represent a candidate element.

    Parameters
    ----------
    node_id: str
        OSM node ID.
    edge_osmid: str
        OSM edge ID.
    obs: tuple
        Coordinate of actual GPS points.
    great_dist: float
        Distance between candidate and actual GPS point.
    coord: tuple
        Candidate coordinate.
    Returns
    -------
    Candidate object
    """
    __slots__ = ["node_id", "edge_osmid", "obs", "great_dist", "coord"]

    def __init__(self, node_id, edge_osmid, obs, great_dist, coord):
        self.node_id = node_id
        self.edge_osmid = edge_osmid
        self.obs = obs
        self.great_dist = great_dist
        self.coord = coord


def prep_candidates(G,interp_dist=10):
    """ Extract candidate points for Hidden-Markov Model map-matching approach.

    Parameters
    ----------
    G: networkx.MultiDiGraph
        Street network graph.
    points: list
        The actual GPS points.
    interp_dist: float, optional, default: 1
        Step to interpolate the graph. The smaller the interp_dist, the greater the precision and the longer
        the computational time.
    closest: bool, optional, default: True
        If true, only the closest point is considered for each edge.
    radius: float, optional, default: 10
        Radius of the search circle.
    Returns
    -------
    G: networkx.MultiDiGraph
        Street network graph.
    results: dict
        Results to be used for the map-matching algorithm.
    """

    if interp_dist:
        G = G.copy()
        if G.graph["geometry"]:
            G = distance.interpolate_graph(G, dist=interp_dist)
        else:
            _ = utils.graph_to_gdfs(G, nodes=False)
            G = distance.interpolate_graph(G, dist=interp_dist)

    geoms = utils.graph_to_gdfs(G, nodes=False).set_index(["u", "v"])[["osmid", "geometry"]]

    uv_xy = [[u, osmid, np.deg2rad(xy)] for uv, geom, osmid in
             zip(geoms.index, geoms.geometry.values, geoms.osmid.values)
             for u, xy in zip(uv, geom.coords[:])]
    
    index, osmid, xy = zip(*uv_xy)
    index = np.array(index)

    y,x = zip(*xy)
    nodes = gpd.GeoSeries.from_xy(x,y,index=osmid)
    nodes_deg = nodes.apply(lambda x: transform(x, np.rad2deg))
    rtree = nodes_deg.sindex
    coords = nodes.get_coordinates()

    ball = BallTree(coords, metric='haversine')

    return index,nodes,nodes_deg,rtree,coords,ball

def get_candidates(G,index,nodes,nodes_deg,rtree,coords,ball, points, closest=True, radius=10):
    """ Extract candidate points for Hidden-Markov Model map-matching approach.

    Parameters
    ----------
    G: networkx.MultiDiGraph
        Street network graph.
    points: list
        The actual GPS points.
    interp_dist: float, optional, default: 1
        Step to interpolate the graph. The smaller the interp_dist, the greater the precision and the longer
        the computational time.
    closest: bool, optional, default: True
        If true, only the closest point is considered for each edge.
    radius: float, optional, default: 10
        Radius of the search circle.
    Returns
    -------
    G: networkx.MultiDiGraph
        Street network graph.
    results: dict
        Results to be used for the map-matching algorithm.
    """
    idxs, dists = ball.query_radius(np.deg2rad(points), radius / distance.EARTH_RADIUS_M, return_distance=True)
    not_found = [[points[i],i] for i, idx_array in enumerate(idxs) if idx_array.size == 0]
    while len(not_found) > 0:
        not_found, nf_i = zip(*not_found)
        radius = radius * 2
        nf_idxs, nf_dists = ball.query_radius(np.deg2rad(not_found), radius / distance.EARTH_RADIUS_M, return_distance=True)
        idxs[list(nf_i)] = nf_idxs
        dists[list(nf_i)] = nf_dists
        not_found = [[points[i],i] for i, idx_array in enumerate(idxs) if idx_array.size == 0]
    
    dists = dists * distance.EARTH_RADIUS_M  # radians to meters

    if closest:
        results = dict()
        for i, (point, idx, dist) in enumerate(zip(points, idxs, dists)):
            df = pd.DataFrame(
                            {"osmid": list(index[idx]),
                               "edge_osmid": list(nodes.index[idx]),
                               "coords": tuple(map(tuple, np.rad2deg(coords.values[idx]))),
                               "dist": dist})

            df = df.loc[df.groupby('edge_osmid')['dist'].idxmin()].reset_index(drop=True)

            results[i] = {
                "observation": point,
                          "osmid": list(df["osmid"]),
                          "edge_osmid": list(df["edge_osmid"]),
                          "candidates": list(df["coords"]),
                          "candidate_type": np.full(len(df["coords"]), False),
                          "dists": list(df["dist"])}

    else:
        results = {i: {
            "observation": point,
                       "osmid": list(index[idx]),
                       "edge_osmid": list(nodes.index[idx]),
                       "candidates": list(map(tuple, np.rad2deg(coords.values[idx]))),
                       "candidate_type": np.full(len(nodes.index[idx]), False),
                       "dists": list(dist)} for i, (point, idx, dist) in enumerate(zip(points, idxs, dists))}

    no_cands = [node_id for node_id, cand in results.items() if not cand["candidates"]]

    if no_cands:
        for cand in no_cands:
            del results[cand]
        print(f"A total of {len(no_cands)} points has no candidates: {*no_cands,}.",end=' ',flush=True)

    dist = max(max(results.values(), key=lambda point: max(point['dists']))['dists'])
    polygon = extract_unique_points(MultiPoint([latlng for point in results.values() for latlng in point['candidates']]))
    # plotting.plot_points(polygon,ax=ax,color='red')

    poly_buff = polygon.convex_hull.buffer(np.rad2deg(max(1000,dist*10)/distance.EARTH_RADIUS_M),cap_style=3)
    # plotting.plot_polygon(polygon.convex_hull,ax=ax,color='green',alpha=0.7)
    # plotting.plot_line(poly_buff.boundary,ax=ax,color='blue',alpha=0.5)
    # plotting.plot_polygon(MultiPoint(nodes_deg.array).convex_hull,ax=ax,color='yellow',alpha=0.3)
    # plt.show()
    if poly_buff.is_valid and poly_buff.area > 0:
        possible_matches_iloc = rtree.query(poly_buff,predicate='intersects')
        possible_matches = nodes_deg.iloc[list(possible_matches_iloc)]

        precise_matches = possible_matches[possible_matches.intersects(poly_buff)]
        to_keep = set(precise_matches.index)
            
    try:
        to_remove = set(index[~nodes_deg.index.isin(to_keep)])
        to_keep = set(index[nodes_deg.index.isin(to_keep)])
    except NameError as err:
        raise ValueError("Found no graph nodes within the requested polygon.") from err

    if len(to_keep) < len(to_remove):
        G_sub = G.subgraph(to_keep).copy()
    else:
        G_sub = G.copy().remove_nodes_from(to_remove)
    print(f'Graph size reduced by {len(to_remove)/len(G.nodes):.0%}.',end=' ',flush=True)

    return G_sub, results


def elab_candidate_results(results, predecessor):
    """ Elaborate results of ``candidate.get_candidates`` method. It selects which candidate best matches the actual
    GPS points.

    Parameters
    ----------
    results: dict
        Output of ``candidate.get_candidates`` method.
    predecessor: dict
        Output of ``mpmatching.viterbi_search`` method.
    Returns
    -------
    results: dict
        Elaborated results.
    """
    results = copy.deepcopy(results)
    for key in list(results.keys())[:-1]:
        win_cand_idx = int(predecessor[str(key + 1)].split("_")[1])
        results[key]["candidate_type"][win_cand_idx] = True

    win_cand_idx = int(predecessor["target"].split("_")[1])
    results[list(results.keys())[-1]]["candidate_type"][win_cand_idx] = True
    return results
