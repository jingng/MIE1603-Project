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

def get_lower_bound(relaxation_sol, nodes, edges):

    new_solution = relaxation_sol

    # Create solution graph
    G = nx.DiGraph()

    for edge in edges:
        i,j = edge

        if relaxation_sol[edge]>0.05:
            G.add_edge(i,j, visited=relaxation_sol[edge])

    # Remove all edges from G where visited is an integer number
    for edge in G.edges:
        G.edges[edge]['visited'] = G.edges[edge]['visited']%1


    for node in G.nodes:

        in_edges = G.in_edges(node) 
        out_edges = G.out_edges(node)

        G.nodes[node]['demand'] = sum(G.edges[edge]['visited'] for edge in in_edges) - sum(G.edges[edge]['visited'] for edge in out_edges)
        G.nodes[node]['supply'] = sum(G.edges[edge]['visited'] for edge in out_edges) - sum(G.edges[edge]['visited'] for edge in in_edges)


    # Check if the supply at all nodes is equal to 0, ie. if set of nodes where supply!=0 is empty, then all nodes[supply]=0
    supply_nodes = {nodes for nodes in G if G.nodes[node]['supply']!=0}
    demand_nodes = {nodes for nodes in G if G.nodes[node]['demand']!=0}

    if len(supply_nodes) == 0:

        for edge in edges:
            new_solution[edge]=int(relaxation_sol[edge])

    if len(demand_nodes) < 10:
        flow = nx.min_cost_flow(G,weight="visited")

    for u,v in G.edges:
        if u in flow and v in flow[u]:
            new_solution[(u,v)]= flow[u][v]
        else:
            new_solution[(u,v)] = 0

    return new_solution


def callback(model, where):
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
      visited_edges = [edge for edge in model._edges if solution_x[edge]>0.05]
      visited_nodes = set()
      for edge in visited_edges:
        i,j = edge
        visited_nodes.update({i,j})

      # Define graph of subtours
      S = nx.Graph()
      for i in visited_nodes:
         for j in visited_nodes:
            if (i,j) in model._edges:
                if i!=j and solution_x[(i,j)] > 0.05:
                    S.add_edge(i,j)

      # Find connected components (subtours)
      components = list(nx.connected_components(S))

      # Add subtour elimination constraints for each disconnected component
      for subtour in components:
        if model._depot in subtour:
          continue

        expr_arc = gp.quicksum(model._x[i,j] for i in subtour for j in subtour if i!=j and (i,j) in model._edges)
        model.cbLazy((expr_arc <= len(subtour) - 1))

    if where ==  GRB.Callback.MIPNODE and GRB.Callback.MIPNODE_STATUS ==  GRB.OPTIMAL:

        node_count = model.cbGet(GRB.Callback.MIPNODE_NODCNT)

        if node_count<=10 or (node_count <= 50 and node_count % 2 == 0) or node_count % 10 ==0:

            relaxation_sol = model.cbGetNodeRel(model._x)
            new_solution = get_lower_bound(relaxation_sol, model._nodes, model._edges)

            # Objective value of current solution
            current_obj_val = model.cbGet(GRB.Callback.MIPSOL_OBJBST)

            # Evaluating objective value at new solution
            objective = model.getObjective()
            new_sol_obj = 0
            for i,j in objective.getVarCoeff():
                new_solution += new_solution[i]*j

            if new_sol_obj > current_obj_val:
                model.cbSetSolution(model._x, new_solution)

def solve_MIP(C, total_length, start_node_no, time_limit):

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

    ## Constraints
    print('adding constraints...')

    # Symmetry Constraint
    symm = m.addConstrs((x.sum(i,'*') == x.sum('*', i) for i in node_dict), name = 'symmetry')

    # Connectivity Constraint (bi-conditional)
    conn1_const = m.addConstrs((x[i,j] >= y[i,j] for i in nodes for j in nodes if i!=j and j in node_dict[i]))
    conn2_const = m.addConstrs((y[i,j]*2 >= x[i,j] for i in nodes for j in nodes if i!=j and j in node_dict[i]))

    # Total Length Constraint
    max_length_const = m.addConstr((x.prod(lengths)<=total_length), name = 'total_length')
    min_length_const = m.addConstr((x.prod(lengths)>=total_length*0.8), name = 'total_length')


    # Start at start node
    start_node = m.addConstrs((y.sum(start_node_no, j)>=1 for j in node_dict[start_node_no]), name='start_node') ## this crashes the model, makes it infeasible

    # Set only one route as profitable

    profit_1 = m.addConstrs((z[i,j] >= y[i,j] for i in nodes for j in nodes if i<=j and (i,j) in edges))
    profit_2 = m.addConstrs((z[i,j] >= y[j,i] for i in nodes for j in nodes if i<=j and (i,j) in edges))
    profit_3 = m.addConstrs((y[i,j] + y[j,i] >= z[i,j] for i in nodes for j in nodes if i<=j and (i,j) in edges))

    print('add objective function...')
    m.setObjective(z.prod(profit)-z.prod(penalty)-x.prod(visited), GRB.MAXIMIZE)

    # Set-up lazy constraints and call variables
    m.Params.LazyConstraints = 1
    m._x = x
    m._y = y
    m._z = z
    m._edges = edges
    m._nodes = nodes
    m._data = []
    m._depot = start_node_no
    m.setParam('TimeLimit', time_limit)

    print('start solving model...')
    m.optimize(callback)

    # m.write('instance-jess.lp')
    # m.optimize()

    # for v in m.getVars():
    #     if v.x > 1e-6:
    #         print(v.varName, v.x)

    print('Total profit: ', m.objVal)

    # Draw final graph
    F = nx.DiGraph()

    for i in nodes:
        for j in nodes:
            if (i,j) in edges:
                if x[(i,j)].x >= 0.05:
                    F.add_edge(i,j, visited=float(x[(i,j)].x))

    # Find total distance
    run_length = 0
    for edge in F.edges:
        run_length += lengths[edge]*F.edges[edge]['visited']

    print('Total length : ' + str(run_length))

    print('Time to solve :' + str(m.Runtime))

    pos = nx.spring_layout(F)
    plt.figure()
    nx.draw(F, pos, with_labels=True, node_color='lightblue', edge_color='gray', arrows=True)
    plt.title("Optimal Run")
    # plt.show()

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

    with open('sol_graph_'+instance + '_' + str(total_length)+ str('_test_better_lb') +'.pkl', "wb") as file:
        pickle.dump(F, file, pickle.HIGHEST_PROTOCOL)

    print("complete!")


# Parameters
total_length_lst = [1000, 5000, 8000, 10000, 15000]
instance_lst = ["instance-jess-min","instance-jess"]#, "instance-jess"]#, "instance-jin.pkl"]

start_node_no = 6813225352 # Start node for instance-jess-min
# start_node_no = 1004361926 # New start node for instance-jess

time_limit=3600

for instance in instance_lst:
  for total_length in total_length_lst:

    log_filepath = "C://Users//Jin//Github//MIE1603-Project//experiment_jess_3600//"
    log_filename = log_filepath + instance + '_' + str(total_length) + str('_test_better_lb')

    print("reading instance file...")

    C  = pd.read_pickle(instance +'.pkl')

    C = ox.truncate.truncate_graph_dist(C,start_node_no,total_length)

    # Run experiments
    print('running experiment with ' + str(total_length) + ' total_length')
    solve_MIP(C, total_length, start_node_no, time_limit)


