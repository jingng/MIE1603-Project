import itertools
import math

import networkx as nx
from shapely.geometry import LineString

from runtrack.graph import distance
from runtrack.matching import candidate

SIGMA_Z = 4.07
BETA = 3

def _emission_prob(dist, sigma=SIGMA_Z):
    """ Compute emission probability of a node

    Parameters
    ----------
    dist: float
        Distance between a real GPS point and a candidate node.
    sigma: float, optional, default: SIGMA_Z
        It is an estimate of the magnitude of the GPS error. See https://www.ismll.uni-hildesheim.de/lehre/semSpatial-10s/script/6.pdf
        for a more detailed description of its calculation.

    Returns
    -------
    ret: float
        Emission probability of a node.
    """
    c = 1 / (sigma * math.sqrt(2 * math.pi))
    return c * math.exp(-(dist / sigma) ** 2)


# A gaussian distribution
def emission_prob(u, sigma=SIGMA_Z):
    """ Compute emission probability of a node

    Parameters
    ----------
    u: runtrack.matching.Candidate
        Node of the graph.
    sigma: float, optional, default: SIGMA_Z
        It is an estimate of the magnitude of the GPS error. See https://www.ismll.uni-hildesheim.de/lehre/semSpatial-10s/script/6.pdf
        for a more detailed description of its calculation.

    Returns
    -------
    ret: float
        Emission probability of a node.
    """
    if u.great_dist:
        c = 1 / (sigma * math.sqrt(2 * math.pi))
        return c * math.exp(-(u.great_dist / sigma) ** 2)
    else:
        return 1

# A empirical distribution
def transition_prob(G, u, v, saved_dists, beta=BETA):
    """ Compute transition probability between node u and v.

    Parameters
    ----------
    G: networkx.MultiDiGraph
        Road network graph.
    u: dict
        Starting node of the graph.
    v: dict
        Target node of the graph.
    beta: float
        This describes the difference between route distances and great circle distances. See https://www.ismll.uni-hildesheim.de/lehre/semSpatial-10s/script/6.pdf
        for a more detailed description of its calculation.

    Returns
    -------
    ret: float
        Transition probability between node u and v.
    """
    c = 1 / beta

    if u.great_dist and v.great_dist:
        s = u if u.node_id < v.node_id else v
        t = u if u.node_id > v.node_id else v
        if s.node_id not in saved_dists or t.node_id not in saved_dists[s.node_id]:
            weight = nx.algorithms.shortest_paths.weighted._weight_function(G, 'length')
            paths = {s.node_id: [s.node_id]}
            sp_dist = nx.algorithms.shortest_paths.weighted._dijkstra(G, s.node_id, weight, paths=paths,target=t.node_id)
            if s.node_id not in saved_dists:
                saved_dists[s.node_id] = {}
            for target, path in reversed(paths.items()):
                if target not in sp_dist:
                    sp_dist[target] = sum(weight(path[i],path[i+1],G[path[i]][path[i+1]]) for i in range(len(path)-1))
                for idx, pred in enumerate(path):
                    if pred not in sp_dist:
                        sp_dist[pred] = path[idx-1] + weight(path[idx-1],pred,G[path[idx-1]][pred])
                    s_path = pred if pred < target else target
                    t_path = pred if pred > target else target
                    if s_path not in saved_dists:
                        saved_dists[s_path] = {}                        
                    if t_path not in saved_dists[s_path]:
                        saved_dists[s_path][t_path] = {'shortest_path': sp_dist[target] - sp_dist[pred]}
                    else:
                        break
            try:
                saved_dists[s.node_id][t.node_id]
            except KeyError as err:
                for pred in sp_dist.keys():
                    if pred < t.node_id:
                        saved_dists[pred][t.node_id] = None
                    else:
                        saved_dists[t.node_id][pred] = None
                        
        try:
            if 'haversine_dist' not in saved_dists[s.node_id][t.node_id]:
                saved_dists[s.node_id][t.node_id]['haversine_dist'] = distance.haversine_dist(*u.coord, *v.coord)
            delta = abs(
                saved_dists[s.node_id][t.node_id]['shortest_path']
                - saved_dists[s.node_id][t.node_id]['haversine_dist'])
        except TypeError as err:
            raise nx.NetworkXNoPath(f"No path to {v.node_id}.") from err
        return c * math.exp(-delta / beta)
    else:
        return 1

def create_trellis(results):
    """ Create a Trellis graph.

    Parameters
    ----------
    results: dict
        Output of ``candidate.get_candidates`` method.
    Returns
    -------
    G: networkx.DiGraph
        A directed acyclic Trellis graph.
    """

    G = nx.DiGraph()

    prev_nodes = ["start"]
    G.add_node("start", candidate=candidate.Candidate("start", None, None, None, None))
    G.add_node("target", candidate=candidate.Candidate("start", None, None, None, None))

    for idx, item in results.items():
        obs = item["observation"]
        nodes, data_node = zip(
            *[(f"{idx}_{j}", {"candidate": candidate.Candidate(node_id, edge_osmid, obs, dist, cand)})
              for j, (node_id, edge_osmid, dist, cand) in
              enumerate(zip(item["osmid"], item["edge_osmid"], item["dists"], item["candidates"]))])
        edges = [(u, v) for u in prev_nodes for v in nodes]

        G.add_nodes_from(zip(nodes, data_node))
        G.add_edges_from(edges)
        prev_nodes = nodes

    edges = [(u, "target") for u in prev_nodes]
    G.add_edges_from(edges)
    return G


def create_path(G, trellis, predecessor):
    """ Create the path that best matches the actual GPS data.

    Parameters
    ----------
    G: networkx.MultiDiGraph
        Road network graph.
    trellis: networkx.DiGraph
        A directed acyclic graph.
    predecessor: dict
        Predecessor for each node.
    Returns
    -------
    path_elab: list
        List of node IDs.
    """
    u, v = list(zip(*[(u, v) for v, u in predecessor.items()][::-1]))
    path = [(u, v) for u, v in zip(u, u[1:])]
    path_elabs = []
    path_elab = []
    isolates = set()

    for (u,v) in path:
        if u.startswith('gap'):
            continue
        if v.startswith('gap'):
            if path_elab:
                path_elabs.append(path_elab)
                path_elab = []
            else:
                isolates.add(trellis.nodes[u]["candidate"].node_id)
            continue

        path_elab.extend(nx.shortest_path(G, trellis.nodes[u]["candidate"].node_id,
                            trellis.nodes[v]["candidate"].node_id,
                            weight='length'))
    path_elabs.append(path_elab)
    for idx,path_elab in enumerate(path_elabs):
        path_elabs[idx] = [k for k, g in itertools.groupby(path_elab)]
        if len(path_elab) == 1:
            isolates.add(path_elab[0])
    path_elabs = [path_elab for path_elab in path_elabs if len(path_elab) > 1]
    if isolates:
        print(f'A total of {len(isolates)} unreachable nodes discarded: {*isolates,}.',end=' ',flush=True)
    else:
        print('All nodes included.',end=' ',flush=True)
    return path_elabs


def create_matched_path(G, trellis, predecessor):
    """ Create the path that best matches the actual GPS points. Route created based on results obtained from ``pmatching_utils.viterbi_search`` and ``mpmatching_utils.create_trellis`` methods.

    Parameters
    ----------
    G: networkx.MultiDiGraph
        Street network graph used to create trellis graph.
    trellis: networkx.DiGraph
        A directed acyclic Trellis graph.
    predecessor: dict
        Predecessor for each node.

    Returns
    -------
    node_ids: list
        List of ids of the nodes that compose the path.
    path_coords: list
        List of nodes' coordinates, in the form of tuple (lat, lon), composing the path.
    """
    node_ids = create_path(G, trellis, predecessor)
    path_coords = [(lat, lng) for lng, lat in LineString([G.nodes[node]["geometry"] for node in node_ids]).coords]
    return node_ids, path_coords


def get_predecessor(target, predecessor):
    """ Reconstruct predecessor dictionary of a decoded trellis DAG.

    Parameters
    ----------
    target: str
        Target node of the trellis DAG.
    predecessor: dict
        Dictionary containing the predecessors of the nodes of a decoded Trellis DAG.
    Returns
    -------
    pred_elab: dict
        Dictionary containing the predecessors of the best nodes of a decoded Trellis DAG.
    """
    pred_elab = {}
    pred = predecessor[target]
    while pred != "start":
        if target.startswith('gap'):
            pred_elab[target] = pred
        else:
            pred_elab[target.split("_")[0]] = pred
        target = pred
        pred = predecessor[target]
    return pred_elab
