import gurobipy as gp
from gurobipy import GRB

import networkx as nx

total_length = 45

n_nodes = 10
nodes_lst = [i for i in range(1, n_nodes+1)]

h_edge_lst = [(i,j) for i in nodes_lst for j in nodes_lst if j!=i] # all connected

# (node i, node j, {"length": 5})
h_edge_lst = [
            (1,2),
            (1,3),
            (1,7),
            (2,4),
            (3,4),
            (2,5),
            (4,6),
            (4,7),
            (7,8),
            (5,9),
            (5,8),
            (6,9),
            (6,10),
            (8,10),
            (9,10)]

h_edge_lst = h_edge_lst + [tuple(reversed(edge)) for edge in h_edge_lst] # adding reversed direction of edges

# Dictionary of connected nodes
node_dict = {i:{j for _, j in h_edge_lst if _ == i } for i in nodes_lst}

g_edge_lst = [
            (1,3),
            (2,4),
            (4,6),
            (4,7),
            (5,9),
            (6,9),
            (1,7),
            (8,10)]

# Define H (runnable paths)
H = nx.DiGraph()
H.add_nodes_from(nodes_lst)
H.add_edges_from(h_edge_lst)
for edge in H.edges:
    H.edges[edge]['length']=5

# list(H.edges(data=True))
# H.edges[1,2]['distance']

# Define G (citystrides)
G = nx.DiGraph()
G.add_nodes_from(nodes_lst)
G.add_edges_from(g_edge_lst)

# 1. All roads have been ran, R = H
R = nx.DiGraph(H)

# Add passes
for edge in R.edges:
    R.edges[edge]['passes']=1

# Define model
m = gp.Model("mp")

edges, lengths = gp.multidict({edge:H.edges[edge]['length'] for edge in H.edges })
profit = {edge: (10 if edge in G.edges else 0) for edge in H.edges}
penalty = {edge: (2 if edge in R.edges else 0) for edge in H.edges}

x = m.addVars(edges, name="traversed", vtype=GRB.INTEGER) # number of times arc (i,j) is traversed
x.update({(j,i):x[i,j] for i,j in edges}) 

y = m.addVars(edges, name="run", vtype=GRB.BINARY) # arc (i,j) is chosen to be ran, {0,1}
y.update({(j,i):y[i,j] for i,j in edges})

## Constraints

# Symmetry Constraint

symm = m.addConstrs((gp.quicksum(x[i,j] for j in node_dict[i] if j != i) == gp.quicksum(x[j,i] for j in node_dict[i] if j != i ) for i in nodes_lst), name='symmetry')

# Connectivity Constraint (Subsets) 

# Connectivity Constraint
conn_const = m.addConstrs((x[i,j] >= y[i,j] for i in nodes_lst for j in nodes_lst if i!=j and j in node_dict[i]))

# Total Length Constraint
total_length_const = m.addConstr((x.prod(lengths)<=total_length), name = 'total_length')

# Start at start node
start_node = m.addConstrs((y.sum(1, j)>=1 for j in node_dict[1]), name='start_node')

m.setObjective(y.prod(profit)-x.prod(penalty), GRB.MAXIMIZE)
# m.write('optimal_run.lp')
m.optimize()

for v in m.getVars():
    if v.x > 1e-6:
        print(v.varName, v.x)

print('Total profit: ', m.objVal)