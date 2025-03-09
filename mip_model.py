import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import matplotlib.pyplot as plt
import osmnx as ox
from collections import defaultdict

import generate_instances
import pickle
import pandas as pd
import time

# Parameters
total_length=100
start_node_no = 8581834456 # Start node for delta, in instance-jess-min and instance-jess graphs
time_limit=300
instance_name = "instance-jess-min"

log_filepath = "C://Users//Jin//Github//MIE1603-Project//logs_cv//"
log_filename = log_filepath + instance_name + '_' + str(total_length)

print("reading instance file...")
C = pd.read_pickle("instance-jess-min.pkl")
# C = pd.read_pickle("instance-jess.pkl")
# C = pd.read_pickle("instance-jin.pkl")

# def solve_MIP(C, total_length, start_node_no, time_limit)

C = nx.DiGraph(C)

# Get visited, profit, and penalty
visited = {edge: C.edges[edge]['visited'] for edge in C.edges}
profit = {edge: C.edges[edge]['profit'] for edge in C.edges}
penalty = {edge: C.edges[edge]['penalty'] for edge in C.edges}

edges, lengths = gp.multidict({edge:C.edges[edge]['length'] for edge in C.edges })
nodes = C.nodes
# node_dict = {i:{j for _, j in edges if _ == i } for i in nodes} # Dictionary of connected nodes, old really slow

node_dict = defaultdict(set)
for i,j in edges:
  node_dict[i].add(j)

# Build model

m = gp.Model("mp")
print('building model...')

# Variables
print('adding variables...')
x = m.addVars(edges, name="traversed", vtype=GRB.INTEGER) # number of times arc (i,j) is traversed
y = m.addVars(edges, name="run", vtype=GRB.BINARY) # arc (i,j) is chosen to be ran, {0,1}

z = m.addVars({(i,j) for i,j in edges if i< j}, name = 'profitable_route', vtype=GRB.BINARY)
# z = m.addVars({(i,j) for i in nodes for j in nodes if i<j}, name='profitable_route', vtype=GRB.BINARY)

## Constraints
print('adding constraints...')

# Symmetry Constraint
symm = m.addConstrs((x.sum(i,'*') == x.sum('*', i) for i in node_dict), name = 'symmetry')

# Connectivity Constraint (bi-conditional)
conn1_const = m.addConstrs((x[i,j] >= y[i,j] for i in nodes for j in nodes if i!=j and j in node_dict[i]))
conn2_const = m.addConstrs((y[i,j]*500 >= x[i,j] for i in nodes for j in nodes if i!=j and j in node_dict[i]))

# Total Length Constraint
total_length_const = m.addConstr((x.prod(lengths)<=total_length), name = 'total_length')

# Start at start node
start_node = m.addConstrs((y.sum(start_node_no, j)>=1 for j in node_dict[start_node_no]), name='start_node') ## this crashes the model, makes it infeasible

# Set only one route as profitable

profit_1 = m.addConstrs((z[i,j] >= y[i,j] for i in nodes for j in nodes if i<=j and (i,j) in edges))
profit_2 = m.addConstrs((z[i,j] >= y[j,i] for i in nodes for j in nodes if i<=j and (i,j) in edges))
profit_3 = m.addConstrs((y[i,j] + y[j,i] >= z[i,j] for i in nodes for j in nodes if i<=j and (i,j) in edges))

print('add objective function...')
m.setObjective(z.prod(profit)-x.prod(penalty)-y.prod(visited), GRB.MAXIMIZE)

def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:

      # Get results
      obj=model.cbGet(GRB.Callback.MIPSOL_OBJBST)
      bound=model.cbGet(GRB.Callback.MIPSOL_OBJBND)
      time=model.cbGet(GRB.Callback.RUNTIME)
      
      model._data.append((time, obj, bound))

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
m._z = z
m._edges = edges
m._N = nodes
m._data = []
m._depot = start_node_no
m.setParam('TimeLimit', time_limit)

print('start solving model...')
m.optimize(subtourelim)

m.write('instance-jess.lp')
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
plt.title("Optimal Run")

# Add last solution data
obj=m.ObjVal
bound=m.ObjBound
time=m.Runtime
m._data.append((time, obj, bound))

# Record results
with open(log_filename + '.csv', "w") as f:
  for sol_data in m._data:
    time, obj, bound = sol_data

    f.write("{}, {}, {}\n".format(time, obj, bound))

# Create final graph solution
print("saving graphs to file", end=" ", flush=True)

with open(instance_name + '_' + str(total_length)+'.pkl', "wb") as file:
    pickle.dump(F, file, pickle.HIGHEST_PROTOCOL)

print("complete!")



# F = nx.MultiDiGraph(F)

# # Plot map 
# fig, ax = ox.plot.plot_graph(
#     F,
#     show=False,
#     close=False,
#     bgcolor="#333333",
#     edge_color="w",
#     edge_linewidth=0.3,
#     node_size=1,
#     node_color='r'
# )

plt.show()


#######################################################################################

