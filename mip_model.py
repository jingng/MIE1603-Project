import gurobipy as gp
from gurobipy import GRB
import networkx as nx

total_length = 45

n_nodes = 9
nodes_lst = [i for i in range(1, n_nodes+1)]

# h_edge_lst = [(i,j) for i in nodes_lst for j in nodes_lst if j!=i] # all connected

# (node i, node j)
h_edge_lst = [
            (1,2),
            (1,6),
            (1,7),
            (2,3),
            (2,7),
            (2,8),
            (3,4),
            (3,8),
            (4,5),
            (4,8),
            (5,6),
            (5,9),
            (6,7),
            (7,8),
            (7,9),
            (8,9)
            ]

h_edge_lst = h_edge_lst + [tuple(reversed(edge)) for edge in h_edge_lst] # adding reversed direction of edges

# Dictionary of connected nodes
node_dict = {i:{j for _, j in h_edge_lst if _ == i } for i in nodes_lst}

g_edge_lst = [
            (1,2),
            (3,4),
            (4,5),
            (5,4),
            (6,1),
            (1,7),
            (7,1),
            (5,9),
            (9,5),
            (5,7),
            (7,9), 
            (8,7)
            ]

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
    R.edges[edge]['visited']=1

# Define model
m = gp.Model("mp")

edges, lengths = gp.multidict({edge:H.edges[edge]['length'] for edge in H.edges })
profit = {edge: (10 if edge in G.edges and R.edges[edge]['visited']==0 else 1) for edge in H.edges}
penalty = {edge: (5 if edge in R.edges else 1) for edge in H.edges}

x = m.addVars(edges, name="traversed", vtype=GRB.INTEGER) # number of times arc (i,j) is traversed
# x.update({(j,i):x[i,j] for i,j in edges}) 

y = m.addVars(edges, name="run", vtype=GRB.BINARY) # arc (i,j) is chosen to be ran, {0,1}
# y.update({(j,i):y[i,j] for i,j in edges})

z = m.addVars({(i,j) for i in nodes_lst for j in nodes_lst if i<j}, name='profitable_route', vtype=GRB.BINARY)

## Constraints

# Symmetry Constraint
symm = m.addConstrs((x.sum(i,'*') == x.sum('*', i) for i in node_dict), name = 'symmetry')

# Connectivity Constraint
conn1_const = m.addConstrs((x[i,j] >= y[i,j] for i in nodes_lst for j in nodes_lst if i!=j and j in node_dict[i]))
conn2_const = m.addConstrs((y[i,j]*10 >= x[i,j] for i in nodes_lst for j in nodes_lst if i!=j and j in node_dict[i]))

# Total Length Constraint
total_length_const = m.addConstr((x.prod(lengths)<=total_length), name = 'total_length')

# Start at start node
start_node = m.addConstrs((y.sum(1, j)>=1 for j in node_dict[1]), name='start_node')

# Set which routes are profitable

profit_1 = m.addConstrs((z[i,j] >= y[i,j] for i in nodes_lst for j in nodes_lst if i<=j and (i,j) in edges))
profit_2 = m.addConstrs((z[i,j] >= y[j,i] for i in nodes_lst for j in nodes_lst if i<=j and (i,j) in edges))
profit_3 = m.addConstrs((y[i,j] + y[j,i] >= z[i,j] for i in nodes_lst for j in nodes_lst if i<=j and (i,j) in edges))

m.setObjective(z.prod(profit)-x.prod(penalty), GRB.MAXIMIZE)

def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:

      solution_x = model.cbGetSolution(model._x)
      solution_y = model.cbGetSolution(model._y)

    #   ## Reverse Arc Constraint
    #   for i in range(1, model._N):
    #     for j in range(1, model._N):

    #       if (i,j) in model._edges and (j,i) in model._edges:
    #         if i!=j and solution_y[(i,j)] > 0.05 and solution_y[(j,i)] > 0.05 :
    #           model.cbLazy(model._y[(i,j)] + model._y[(j,i)] <= 1)

    #   Testing on one set of edges
    #   if solution_x[(1,7)] > 0.05 and solution_x[(7,1)] > 0.05 :
    #       model.cbLazy(model._y[(1,7)] + model._y[(7,1)] <= 1)
    #   if solution_x[(5,9)] > 0.05 and solution_x[(9,5)] > 0.05 :
    #       model.cbLazy(model._y[(5,9)] + model._y[(9,5)] <= 1)
    #   if solution_x[(1,2)] > 0.05 and solution_x[(1,2)] > 0.05 :
    #       model.cbLazy(model._y[(1,2)] + model._y[(2,1)] <= 1)

    ## Subtour elimination
    # Identify visited nodes
      visited_edges = [edge for edge in model._edges if solution_x[edge]>0.05]
      visited_nodes = set()
      for edge in visited_edges:
        i,j = edge
        visited_nodes.update({i,j})

      # Define graph of subtours
      S = nx.Graph()
      for i in visited_nodes:
         for j in visited_nodes:
            if (i,j) in edges:
                if i!=j and solution_x[(i,j)] > 0.05:
                    S.add_edge(i,j)

      # Find connected components (subtours)
      components = list(nx.connected_components(S))

      # Add subtour elimination constraints for each disconnected component
      for subtour in components:
        if 0 in subtour:
          continue

        expr_arc = gp.quicksum(model._y[i,j] for i in subtour for j in subtour if i!=j and (i,j) in model._edges)

        model.cbLazy((expr_arc <= len(subtour) - 1))


# Set-up lazy constraints and call variables
m.Params.LazyConstraints =1
m._x = x
m._y = y
m._edges = edges
m._N = n_nodes
m.optimize(subtourelim)

m.write('optimal_run.lp')
m.optimize()

for v in m.getVars():
    if v.x > 1e-6:
        print(v.varName, v.x)

print('Total profit: ', m.objVal)

