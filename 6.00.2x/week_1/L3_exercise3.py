from nodes_graphs import *

# My graph object
# Create nodes
A = Node("A")
B = Node("B")
C = Node("C")
D = Node("D")

# Create edges
AB = Edge(A, B)
AD = Edge(A, D)
BC = Edge(B, C)
BD = Edge(B, D)
CB = Edge(C, B)

# Create graph
mygraph = Digraph()

# Add nodes to mygraph
nodelist = [A, B, C, D]
for node in nodelist:
    mygraph.addNode(node)

# Add edges to mygraph
edgelist = [AB, AD, BC, BD, CB]
for edge in edgelist:
    mygraph.addEdge(edge)



# g object (from the exercise)
# Generate and add nodes

nodes = []
nodes.append(Node("ABC")) # nodes[0]
nodes.append(Node("ACB")) # nodes[1]
nodes.append(Node("BAC")) # nodes[2]
nodes.append(Node("BCA")) # nodes[3]
nodes.append(Node("CAB")) # nodes[4]
nodes.append(Node("CBA")) # nodes[5]

g = Graph()
for n in nodes:
    g.addNode(n)
    
# Generate and add edges
for n in nodes:
    name = n.getName()
    name1 = name[0] + name[2] + name[1]
    name2 = name[1] + name[0] + name[2]
    for n_dest in nodes:
        if n_dest.getName() == name1:
            g.edges[n].append(n_dest)
        if n_dest.getName() == name2:
            g.edges[n].append(n_dest)

# Define DFS algorithm
def DFS(graph, start, end, path, shortest):
    count = 0
    path = path + [start]
    print("Current DFS path:", [str(node) for node in path])
    if start == end:
        print("A path found:", [str(node) for node in path])
        return path
    for node in graph.childrenOf(start):
        count += 1
        print(f"    Iteration {count} for child node {node}")
        if node not in path:
            if shortest == None or len(path) < len(shortest):
                newPath = DFS(graph, node, end, path, shortest)
            if newPath != None:
                shortest = newPath
                print("        A new short path found:", [str(node) for node in shortest])
        else:
            print("        Node already in path.")
    return shortest

# Define wrapper function for DFS
def shortestPath(graph, start, end):
    return DFS(graph, start, end, [], None)

shortest_path = shortestPath(g, nodes[0], nodes[2])
# shortest_path = shortestPath(mygraph, A, D)

print([str(node) for node in shortest_path])

# print(mygraph)