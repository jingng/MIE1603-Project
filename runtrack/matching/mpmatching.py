import math
from collections import deque
import networkx as nx

from . import mpmatching_utils
from . import candidate


def viterbi_search(G, trellis, saved_dists, start="start", target="target", beta=mpmatching_utils.BETA,
                   sigma=mpmatching_utils.SIGMA_Z):
    """ Function to compute viterbi search and perform Hidden-Markov Model map-matching.

    Parameters
    ----------
    G: networkx.MultiDiGraph
        Street network graph.
    trellis:
    start: str, optional, default: "start"
        Starting node.
    target: str, optional, default: "target"
        Target node.
    beta: float
        This describes the difference between route distances and great circle distances. See https://www.ismll.uni-hildesheim.de/lehre/semSpatial-10s/script/6.pdf
        for a more detailed description of its calculation.

    sigma: float
        It is an estimate of the magnitude of the GPS error. See https://www.ismll.uni-hildesheim.de/lehre/semSpatial-10s/script/6.pdf
        for a more detailed description of its calculation.

    Returns
    -------
    joint_prob: dict
        Joint probability for each node.
    predecessor: dict
        Predecessor for each node.

    Notes
    -----
    See https://www.ismll.uni-hildesheim.de/lehre/semSpatial-10s/script/6.pdf for a more detailed description of this
    method.
    """

    # Initialize joint probability for each node
    joint_prob = {}
    for u_name in trellis.nodes():
        joint_prob[u_name] = -float('inf')
    predecessor = {}
    predecessor_val = {}

    queue = deque()

    queue.append(start)
    joint_prob[start] = math.log10(mpmatching_utils.emission_prob(trellis.nodes[start]["candidate"], sigma))
    predecessor[start] = None
    g = 0
    cnt = 0
    num_nodes = len(trellis)
    logs = [False] * 9

    while queue:
        # Extract node u
        u_name = queue.popleft()
        u = trellis.nodes[u_name]["candidate"]
        bfs = None
        succ_set = set()

        if u_name == target:
            break
        for v_name in list(trellis.successors(u_name)):
            v = trellis.nodes[v_name]["candidate"]

            # if bfs is not None and v_name.split('_')[0] not in ('gap','target') and v.node_id not in succ_set :
            #     for _,successors in bfs:
            #         succ_set.update(successors)
            #         if v.node_id in successors:
            #             break
            #     else:
            #         if not queue:
            #             if u_name.startswith('gap') or u_name == 'start':
            #                 for v_succ in trellis.successors(v_name):
            #                     trellis.add_edge(u_name,v_succ)
            #                 queue.append(u_name)
            #                 continue
            #             trellis.add_node(f'gap_{g}', candidate=candidate.Candidate("start", None, None, None, None))
            #             joint_prob[f'gap_{g}'] = -float('inf')
            #             for v_sib in trellis.successors(u_name):
            #                 trellis.add_edge(f'gap_{g}',v_sib)
            #             for u_sib in trellis.successors(predecessor[u_name]):
            #                 trellis.add_edge(u_sib,f'gap_{g}')
            #                 if u_sib in predecessor.keys():
            #                     queue.append(u_sib)
            #             g = g+1

            try:
                new_prob = joint_prob[u_name] + math.log10(mpmatching_utils.emission_prob(v, sigma)) \
                           + math.log10(mpmatching_utils.transition_prob(G, u, v, saved_dists, beta))

                if joint_prob[v_name] < new_prob:
                    joint_prob[v_name] = new_prob
                    if v_name not in predecessor:
                        predecessor[v_name] = u_name
                        predecessor_val[v_name] = new_prob
                    elif v_name in predecessor and predecessor_val[v_name] < new_prob:
                        predecessor[v_name] = u_name
                        predecessor_val[v_name] = new_prob
                if v_name not in queue:
                    queue.append(v_name)

            except nx.NetworkXNoPath:
                # if bfs is None and list(trellis.successors(u_name)).index(v_name) < len(list(trellis.successors(u_name)))-2:
                #     bfs = nx.bfs_successors(G,u.node_id)
                pass
            except Exception as error:
                pass
                # print('--',error)
        
        if not queue:
            if u_name.startswith('gap') or u_name == 'start':
                for v_succ in trellis.successors(v_name):
                    trellis.add_edge(u_name,v_succ)
                queue.append(u_name)
                continue
            trellis.add_node(f'gap_{g}', candidate=candidate.Candidate("start", None, None, None, None))
            joint_prob[f'gap_{g}'] = -float('inf')
            for v_sib in trellis.successors(u_name):
                trellis.add_edge(f'gap_{g}',v_sib)
            for u_sib in trellis.successors(predecessor[u_name]):
                trellis.add_edge(u_sib,f'gap_{g}')
                if u_sib in predecessor.keys():
                    queue.append(u_sib)
                    cnt = cnt-1
            g = g+1
        
        if not u_name.startswith('gap'):
            cnt = cnt+1
        if False in logs and cnt/num_nodes > 1/10*(logs.index(False)+1):
            logs[logs.index(False)] = True
            print(".", end=" ", flush=True)

    if g > 0:
        print(f'Gaps in route: {g}.',end=' ',flush=True)
    else:
        print('No gaps in route.',end=' ',flush=True)
    predecessor = mpmatching_utils.get_predecessor("target", predecessor)

    return joint_prob[target], predecessor
