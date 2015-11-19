from collections import defaultdict
import Queue as queue

class Vertex:
    WHITE = "white"
    GREY = "grey"
    BLACK = "black"
    def __init__(self, name):
        self.name = name
        self.c = Vertex.WHITE #color
        self.p = None #predecessor
        self.d = float("inf") #distance from source

    def info(self):
        return "{}, {}".format(self.name, self.c)

class Graph:
    def __init__(self):
        self.vertices = []
        self.adj = defaultdict(set) #adj-list for all vertices

    def addEdge(self, u, v, undirected=False):
        self.adj[u].add(v)
        if undirected:
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


def BFS(graph, s):
    # initialize graph for search
    for v in graph.vertices:
        v.c = Vertex.WHITE
        v.d = float("inf")
        v.p = None
    graph.info()

    # begin search
    s.c = Vertex.GREY
    s.d = 0
    s.p = None
    q = queue.Queue()
    q.put(s)
    while not q.empty():
        u = q.get()
        print("Searching edges off {}".format(u.name))
        for v in graph.adj[u]:
            print("    {}".format(v.name))
            if v.c == Vertex.WHITE:
                v.c = Vertex.GREY
                v.d = u.d + 1
                v.p = u
                q.put(v)
        u.c = Vertex.BLACK
        graph.info()

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

main()
