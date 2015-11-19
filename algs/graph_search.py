from collections import defaultdict
import Queue as queue

DEBUG = False

class Vertex:
    WHITE = "white"
    GREY = "grey"
    BLACK = "black"
    def __init__(self, name):
        self.name = name
        self.c = Vertex.WHITE # color to be used in BFS
        self.p = None # predecessor in graph
        self.d = float("inf") # distance from source vertex
        self.f = 0 # timestamp for when finished search, used in DFS

    def info(self, recurse=True):
        info = "{}, {}, {}, {}".format(self.name, self.c, self.d, self.f)
        if self.p is not None and recurse:
            info += ", predecessor: {}".format(self.p.info(False))
        return info

# graph implementation using adj-lists
class Graph:
    def __init__(self):
        self.vertices = [] # ordered list of unique vertices
        self.adj = defaultdict(set)

    def addEdge(self, u, v, undirected=False):
        # add edge to adj list of first vertex
        self.adj[u].add(v)
        if undirected:
            # if undirected, add edge to adj list of second vertex
            self.adj[v].add(u)
        if u not in self.vertices:
            self.vertices.append(u)
        if v not in self.vertices:
            self.vertices.append(v)

    def info(self):
        print("Graph, {} vertices, {} adjMatrix".\
              format(len(self.vertices), len(self.adj)))
        for v in self.vertices:
            print("- {}".format(v.info()))
            for u in self.adj[v]:
                print("  {}".format(u.info()))


# CLRS pg. 595
# breadth-first search where all edges of a vertex are explored at
# even depth levels before going deeper into any one path
def BFS(graph, s):
    # initialize graph for search
    for v in graph.vertices:
        v.c = Vertex.WHITE
        v.d = float("inf")
        v.p = None
    if DEBUG: graph.info()
    # begin search
    s.c = Vertex.GREY
    s.d = 0
    s.p = None
    q = queue.Queue()
    q.put(s)
    while not q.empty():
        u = q.get()
        if DEBUG: print("Searching edges off {}".format(u.name))
        for v in graph.adj[u]:
            if DEBUG: print("    {}".format(v.name))
            if v.c == Vertex.WHITE:
                v.c = Vertex.GREY
                v.d = u.d + 1
                v.p = u
                q.put(v)
        u.c = Vertex.BLACK
        if DEBUG: graph.info()

# CLRS pg. 601
# prints out the vertices on a shortest path from s to v
# assuming that BFS has already computed a breadth-first tree
def PRINT_PATH(graph, s, v):
    if v == s: #if v in graph.adj[s]:
        print(s.info())
    elif v.p == None:
        print("No path from {} to {} exists".format(s.name, v.name))
    else:
        PRINT_PATH(graph, s, v.p)
        print(v.info())

time = 0
# CLRS pg. 604
# depth-first search of graph where path of edge is explored
# to its deepest extent before the next edge path is started
# running time of THETA(V + E)
def DFS(graph):
    global time
    # initialize graph for searching
    for u in graph.vertices:
        u.c = Vertex.WHITE
        u.p = None
    # begin search
    time = 0
    for u in graph.vertices:
        if u.c == Vertex.WHITE:
            DFS_VISIT(graph, u) # called in total |V| times, one per vertex

def DFS_VISIT(graph, u):
    global time
    print("VISITING {}".format(u.name))
    time += 1 # white vertex u has just been discovered
    u.d = time
    u.c = Vertex.GREY
    for v in graph.adj[u]: # explore edges of u, in total |E| times for total sum of all edges in a graph
        if v.c == Vertex.WHITE:
            v.p = u
            DFS_VISIT(graph, v) # immediately visit edges of found edge
    u.c = Vertex.BLACK # blacken u; it is finished
    time += 1
    u.f = time # integer between 1 and 2|V|, since there is only one discovery event and one finishing event for each |V| vertices

def main():
    r = Vertex('r')
    s = Vertex('s')
    t = Vertex('t')
    u = Vertex('u')
    v = Vertex('v')
    w = Vertex('w')
    x = Vertex('x')
    y = Vertex('y')
    z = Vertex('z')

    # graph based on CLRS Fig. 22.3, pg. 596
    graph = Graph()
    graph.addEdge(r, s, True)
    graph.addEdge(r, v, True)
    graph.addEdge(s, w, True)
    graph.addEdge(w, t, True)
    graph.addEdge(w, x, True)
    graph.addEdge(t, x, True)
    graph.addEdge(t, u, True)
    graph.addEdge(u, x, True)
    graph.addEdge(u, y, True)
    graph.addEdge(x, y, True)

    print("===========")
    print("Breadth-First Search")
    print("===========")
    BFS(graph, s)
    PRINT_PATH(graph, s, v)

    # graph based on CLRS Fig. 22.5, pg. 607
    graph2 = Graph()
    graph2.addEdge(s, z)
    graph2.addEdge(s, w)
    graph2.addEdge(z, y)
    graph2.addEdge(z, w)
    graph2.addEdge(t, u)
    graph2.addEdge(t, v)
    graph2.addEdge(u, t)
    graph2.addEdge(u, v)
    graph2.addEdge(v, s)
    graph2.addEdge(v, w)
    graph2.addEdge(w, x)
    graph2.addEdge(y, x)
    graph2.addEdge(x, z)

    print("===========")
    print("Depth-First Search")
    print("===========")
    DFS(graph2)
    PRINT_PATH(graph2, s, v)
    graph2.info()

main()
