# reference, http://www.cse.unt.edu/~tarau/teaching/AnAlgo/Ford%E2%80%93Fulkerson%20algorithm.pdf
import os
import sys

DEBUG = False

class Edge:
    def __init__(self, u, v, capacity):
        self.u = u
        self.v = v
        self.c = capacity

    def __str__(self):
        return "{0}->{1}, c={2}".format(self.u, self.v, self.c)

class NetworkFlowGraph:
    def __init__(self):
        self.adj = {}
        self.flow = {}
        self.residual_flow = {}

    def add_vertex(self, vertex):
        if vertex not in self.adj:
            self.adj[vertex] = []
            self.flow[vertex] = {}
            self.residual_flow[vertex] = {}

    def add_edge(self, u, v, capacity):
        self.add_vertex(u)
        self.add_vertex(v)
        self.adj[u].append(Edge(u, v, capacity))
        # update flow network
        self.flow[u][v] = 0
        self.residual_flow[u][v] = capacity

    def add_source(self, vertex):
        dummy_source = 'source'
        self.add_vertex(dummy_source)
        self.add_vertex(vertex)
        self.add_edge(dummy_source, vertex, float("inf"))

    def add_sink(self, vertex):
        dummy_sink = 'sink'
        self.add_vertex(dummy_sink)
        self.add_vertex(vertex)
        self.add_edge(vertex, dummy_sink, float("inf"))

    def reset_flow(self):
        for u in self.adj:
            for e in self.adj[u]:
                self.flow[e.u][e.v] = 0
                self.residual_flow[e.u][e.v] = e.c

    def add_flow(self, edge, flow):
        self.flow[edge.u][edge.v] += flow
        self.residual_flow[edge.u][edge.v] -= flow

    # get augmenting path via depth first search
    def get_augmenting_path(self,
            source='source', sink='sink', path=[], path_set=set()):
        if source == sink:
            return path
        for e in self.adj[source]:
            residual = self.residual_flow[e.u][e.v]
            if residual > 0 and e not in path_set:
                path_set.add(e)
                if DEBUG: print "adding edge to aug path " + str(e) + \
                        " with remaining flow " + str(residual)
                aug_path = self.get_augmenting_path(
                        e.v, sink, path + [e], path_set)
                if aug_path is not None:
                    return aug_path

    # get augmenting path via breadth first search
    def get_shortest_augmenting_path(self):
        # TODO: implement for EK algorithm!
        pass

    def get_total_flow(self, vertex='source'):
        total_flow = 0
        for v in self.flow[vertex]:
            total_flow += self.flow[vertex][v]
        return total_flow

    def get_min_capacity(self, path):
        min_cap = float("inf")
        for e in path:
            if self.residual_flow[e.u][e.v] < min_cap:
                min_cap = self.residual_flow[e.u][e.v]
        return min_cap

def max_flow_ford_fulkerson(multigraph):
    multigraph.reset_flow()
    path = multigraph.get_augmenting_path()
    while path is not None:
        min_c = multigraph.get_min_capacity(path)
        if DEBUG: print "adding flow " + str(min_c)
        for e in path:
            multigraph.add_flow(e, min_c)
        path = multigraph.get_augmenting_path(path_set=set())
    return multigraph.get_total_flow()

def max_flow_edmonds_karp(multigraph):
    multigraph.reset_flow()
    path = multigraph.get_augmenting_path()
    while path is not None:
        min_c = multigraph.get_min_capacity(path)
        if DEBUG: print "adding flow " + str(min_c)
        for e in path:
            multigraph.add_flow(e, min_c)
        path = multigraph.get_shortest_augmenting_path()
    return multigraph.get_total_flow()

def main():
    network = NetworkFlowGraph()

    network.add_source('x1')
    network.add_source('x2')
    network.add_sink('y1')
    network.add_sink('y2')
    network.add_sink('y3')

    network.add_edge('x1', 'a', 7)
    network.add_edge('x1', 'b', 18)
    network.add_edge('x1', 'c', 5)
    network.add_edge('x2', 'a', 8)
    network.add_edge('x2', 'c', 2)
    network.add_edge('x2', 'd', 6)
    network.add_edge('a', 'y1', 4)
    network.add_edge('a', 'e', 3)
    network.add_edge('a', 'c', 22)
    network.add_edge('b', 'y1', 2)
    network.add_edge('b', 'y2', 7)
    network.add_edge('b', 'a', 19)
    network.add_edge('c', 'e', 24)
    network.add_edge('c', 'y3', 4)
    network.add_edge('d', 'y2', 12)
    network.add_edge('d', 'y3', 15)
    network.add_edge('d', 'b', 13)
    network.add_edge('d', 'c', 24)
    network.add_edge('e', 'b', 7)
    network.add_edge('e', 'd', 16)

    max_flow = max_flow_ford_fulkerson(network)
    print "Max flow from FF: {0}".format(max_flow)
    #max_flow = max_flow_edmonds_karp(network)
    #print "Max flow from EK: {0}".format(max_flow)

main()

