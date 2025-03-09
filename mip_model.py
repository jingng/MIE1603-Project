import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import matplotlib.pyplot as plt
import osmnx as ox

import generate_instances
import pickle
import pandas as pd
import time

C = pd.read_pickle("instance-jess-min.pkl")

# Get visited, profit, and penalty
visited = {edge: C.edges[edge]['visited'] for edge in C.edges}
profit = {edge: C.edges[edge]['profit'] for edge in C.edges}
penalty = {edge: C.edges[edge]['penalty'] for edge in C.edges}

# Plot map 
fig, ax = ox.plot.plot_graph(
    C,
    show=False,
    close=False,
    bgcolor="#333333",
    edge_color="w",
    edge_linewidth=0.3,
    node_size=1,
    node_color='r'
)
# plt.show()

C = nx.DiGraph(C)

# Parameters

total_length=100

edges, lengths = gp.multidict({edge:C.edges[edge]['length'] for edge in C.edges })
nodes = C.nodes

# Dictionary of connected nodes
# node_dict = {i:{(j,k) for _, j, k in edges if _ == i } for i in nodes}
node_dict = {i:{j for _, j in edges if _ == i } for i in nodes}

# Get visited, profit, and penalty
visited = {edge: C.edges[edge]['visited'] for edge in C.edges} # historic penalty
profit = {edge: C.edges[edge]['profit'] for edge in C.edges}
penalty = {edge: C.edges[edge]['penalty'] for edge in C.edges}

start_node_no = 8581834456

# Checking profit
profit = {edge: C.edges[edge]['profit'] for edge in C.edges if C.edges[edge]['profit'] >0}


m = gp.Model("mp")

x = m.addVars(edges, name="traversed", vtype=GRB.INTEGER) # number of times arc (i,j) is traversed
y = m.addVars(edges, name="run", vtype=GRB.BINARY) # arc (i,j) is chosen to be ran, {0,1}

z = m.addVars({(i,j) for i,j in edges if i< j}, name = 'profitable_route', vtype=GRB.BINARY)
# z = m.addVars({(i,j) for i in nodes for j in nodes if i<j}, name='profitable_route', vtype=GRB.BINARY)

## Constraints

# Symmetry Constraint
symm = m.addConstrs((x.sum(i,'*') == x.sum('*', i) for i in node_dict), name = 'symmetry')

# Connectivity Constraint (bi-conditional)
conn1_const = m.addConstrs((x[i,j] >= y[i,j] for i in nodes for j in nodes if i!=j and j in node_dict[i]))
conn2_const = m.addConstrs((y[i,j]*10 >= x[i,j] for i in nodes for j in nodes if i!=j and j in node_dict[i]))

# Total Length Constraint
total_length_const = m.addConstr((x.prod(lengths)<=total_length), name = 'total_length')

# Start at start node
# start_node = m.addConstrs((y.sum(start_node_no, j)>=1 for j in node_dict[start_node_no]), name='start_node') ## this crashes the model, makes it infeasible

# Set only one route as profitable

profit_1 = m.addConstrs((z[i,j] >= y[i,j] for i in nodes for j in nodes if i<=j and (i,j) in edges))
profit_2 = m.addConstrs((z[i,j] >= y[j,i] for i in nodes for j in nodes if i<=j and (i,j) in edges))
profit_3 = m.addConstrs((y[i,j] + y[j,i] >= z[i,j] for i in nodes for j in nodes if i<=j and (i,j) in edges))

m.setObjective(z.prod(profit)-x.prod(penalty)-y.prod(visited), GRB.MAXIMIZE)

def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:

      solution_x = model.cbGetSolution(model._x)
      solution_y = model.cbGetSolution(model._y)

    ## Subtour elimination
    # Identify visited nodes
      visited_edges = [edge for edge in model._edges if solution_y[edge]>0.05]
      visited_nodes = set()
      for edge in visited_edges:
        i,j = edge
        visited_nodes.update({i,j})

      # Define graph of subtours
      S = nx.Graph()
      for i in visited_nodes:
         for j in visited_nodes:
            if (i,j) in edges:
                if i!=j and solution_y[(i,j)] > 0.05:
                    S.add_edge(i,j)

      # Find connected components (subtours)
      components = list(nx.connected_components(S))

      # Add subtour elimination constraints for each disconnected component
      for subtour in components:
        if model._depot in subtour:
          continue

        expr_arc = gp.quicksum(model._y[i,j] for i in subtour for j in subtour if i!=j and (i,j) in model._edges)
        model.cbLazy((expr_arc <= len(subtour) - 1))

# Set-up lazy constraints and call variables
m.Params.LazyConstraints = 1
m._x = x
m._y = y
m._edges = edges
m._N = nodes
m._depot = start_node_no
m.optimize(subtourelim)

# m.write('optimal_run.lp')
# m.optimize()

for v in m.getVars():
    if v.x > 1e-6:
        print(v.varName, v.x)

print('Total profit: ', m.objVal)

# Draw final graph
F = nx.DiGraph()

for i in nodes:
   for j in nodes:
      if (i,j) in edges:
        if x[(i,j)].x >= 0.05:
            F.add_edge(i,j)

pos = nx.spring_layout(F)
plt.figure()
nx.draw(F, pos, with_labels=True, node_color='lightblue', edge_color='gray', arrows=True)
# nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')

plt.title("Optimal Run")
plt.show()


#######################################################################################

