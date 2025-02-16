import networkx as nx

# H: Runnable paths
# G: City Strides Graph, G \subset H
# R: H with number of passes

nodes_lst = [i for i in range(10)]

# (node i, node j, distance)
h_edge_lst = [
            (1,2),
            (1,3),
            (2,4),
            (3,4),
            (2,5),
            (4,6),
            (4,7),
            (7,8),
            (5,9),
            (5,8),
            (6,9),
            (1,7),
            (8,9),
            (6,10),
            (8,10),
            (9,10)]

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

# Define G (citystrides)
G = nx.DiGraph()
G.add_nodes_from(nodes_lst)
G.add_edges_from(g_edge_lst)

# 1. All roads have been ran, R = H
R = nx.DiGraph(H)

# Add passes
for edge in R.edges:
    R.edges[edge]['passes']=1

# 2. No Roads have been ran, R is empty
R = nx.DiGraph()

# 3. Only H/G ran, R = H/G
R = nx.DiGraph(H)
for edge in G.edges:
    R.remove_edge(edge[0], edge[1])

for edge in R.edges:
    R.edges[edge]['passes']=1

# 4. Only G ran, R = G
R = nx.DiGraph(G)

for edge in R.edges:
    R.edges[edge]['passes']=1

# Realistic ones

# 5. Ran 1->7->8->10
r1_edge_lst = [
    (1,7),
    (7,8),
    (8,10)
]

R = nx.DiGraph()
R.add_nodes_from(nodes_lst)
R.add_edges_from(r1_edge_lst)

for edge in R.edges:
    R.edges[edge]['passes']=1

# 6. Ran 1->2->5->9->10->8->7->1
r1_edge_lst = [
    (1,2),
    (2,5),
    (5,9),
    (9,10),
    (7,8),
    (1,7)
]

R = nx.DiGraph()
R.add_nodes_from(nodes_lst)
R.add_edges_from(r1_edge_lst)

for edge in R.edges:
    R.edges[edge]['passes']=1

'''
# Incorporate nodes from one graph into another
H = nx.path_graph(10)
G.add_nodes_from(H)

list(G.nodes) # List all nodes
list(G.edges) # List all edges
list(G.adj[2]) # List adjacent nodes
G.degree[2] # number of edges incident to 1
'''
