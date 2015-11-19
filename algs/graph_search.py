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

    def info(self, recurse=True):
        info = "{}, {}".format(self.name, self.c)
        if self.p is not None and recurse:
            info += ", predecessor: {}".format(self.p.info(False))
        return info

# graph implementation using adj-lists
class Graph:
    def __init__(self):
        self.vertices = set() # ordered list of vertices
        self.adj = defaultdict(set) # adj-list for all vertex edges

    def addEdge(self, u, v, undirected=False):
        self.adj[u].add(v)
        if undirected:
            self.adj[v].add(u)
        self.vertices.add(u)
        self.vertices.add(v)

    def info(self):
        print("Graph, {} vertices, {} adjMatrix".\
              format(len(self.vertices), len(self.adj)))
        for v in self.vertices:
            print("- {}".format(v.info()))
            for u in self.adj[v]:
                print("  {}".format(u.info()))


# CLRS pg. 595
# assumes that the input graph is represented using adj-lists
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
def print_path(graph, s, v):
    if v == s: #if v in graph.adj[s]:
        print(s.info())
    elif v.p == None:
        print("No path from {} to {} exists".format(s.name, v.name))
    else:
        print_path(graph, s, v.p)
        print(v.info())

def main():
    r = Vertex('r')
    s = Vertex('s')
    t = Vertex('t')
    u = Vertex('u')
    v = Vertex('v')
    w = Vertex('w')
    x = Vertex('x')
    y = Vertex('y')

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

    BFS(graph, s)
    print_path(graph, s, v)

main()
