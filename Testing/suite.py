import sys
import numpy as np
from enum import IntEnum
import heapq as heap

class GraphWeight(IntEnum):
    UNWEIGHTED = 0
    WEIGHTED_INT = 1
    WEIGHTED_REAL = 2

class GraphDirectedness(IntEnum):
    UNDIRECTED = 0
    DIRECTED = 1

class Graph(object):
    def __init__(self, n : int, m : int, edges : tuple, directedness : GraphDirectedness, weight : GraphWeight):

        assert m == len(edges)

        self.n = n
        self.m = m
        self.directedness = directedness
        self.weight = weight

        if directedness == GraphDirectedness.UNDIRECTED:
            self.m *= 2
            anti_edges = []
            for e in edges:
                u = e[0]
                v = e[1]
                if weight != GraphWeight.UNWEIGHTED:
                    anti_edges.append((v,u,e[2]))
                else:
                    anti_edges.append((v,u))
            edges += anti_edges

        self.ir = np.zeros(self.n+1, dtype=int)
        self.jc = np.zeros(self.m, dtype=int)

        if self.weight == GraphWeight.UNWEIGHTED:
            self.num = None
        elif self.weight == GraphWeight.WEIGHTED_INT:
            self.num = np.zeros(self.m, dtype=int)
        elif self.weight == GraphWeight.WEIGHTED_REAL:
            self.num = np.zeros(self.m, dtype=float)

        rowcnts = np.zeros(self.n, dtype=int)

        for k in range(self.m):
            u = edges[k][0]
            rowcnts[u] += 1

        nzcnt = 0

        for k in range(self.n):
            self.ir[k] = nzcnt
            nzcnt += rowcnts[k]
            rowcnts[k] = self.ir[k]

        self.ir[self.n] = nzcnt

        for k in range(self.m):
            u, v = edges[k][0:2]
            p = rowcnts[u]
            rowcnts[u] += 1
            self.jc[p] = v
            if self.weight == GraphWeight.WEIGHTED_INT:
                self.num[p] = int(edges[k][2])
            elif self.weight == GraphWeight.WEIGHTED_REAL:
                self.num[p] = float(edges[k][2])

    def is_weighted(self):
        return self.weight != GraphWeight.UNWEIGHTED

    def is_directed(self):
        return self.directedness == GraphDirectedness.DIRECTED

    def is_undirected(self):
        return self.directedness == GraphDirectedness.UNDIRECTED

    def num_vertices(self):
        return self.n

    def num_edges(self):
        return self.m if self.is_directed() else int(self.m/2)

    def write(self, filename):
        with open(filename, "w") as f:
            f.write("{} {} {} {}\n".format(self.num_vertices(), self.num_edges(), self.directedness, self.weight))
            for i in range(self.n):
                for p in range(self.ir[i], self.ir[i+1]):
                    j = self.jc[p]
                    if self.is_undirected() and i > j:
                        continue
                    if self.is_weighted():
                        v = int(self.num[p]) if self.weight == GraphWeight.WEIGHTED_INT else float(self.num[p])
                        f.write("{} {} {}\n".format(i, j, v))
                    else:
                        f.write("{} {}\n".format(i, j))


    @classmethod
    def read(cls, filename):
        with open(filename, "r") as f:
            n, m, directedness, weight = (int(v) for v in next(f).rstrip().split())
            edges = []
            for line in f.readlines():
                items = line.rstrip().split()
                u = int(items[0])
                v = int(items[1])
                if len(items) == 3:
                    w = float(items[2])
                    edges.append((u,v,w))
                else:
                    edges.append((u,v))

        return cls(n, m, edges, directedness, weight)

    def neighbors(self, u, add_weights=False):
        if add_weights:
            if self.is_weighted():
                for p in range(self.ir[u], self.ir[u+1]):
                    yield self.jc[p], self.num[p]
            else:
                for p in range(self.ir[u], self.ir[u+1]):
                    yield self.jc[p], 1
        else:
            for p in range(self.ir[u], self.ir[u+1]):
                yield self.jc[p]

    def edges(self, add_weights=False):
        n = self.num_vertices()
        if add_weights:
            if self.is_weighted():
                for i in range(n):
                    for p in range(self.ir[i], self.ir[i+1]):
                        yield (i, self.jc[p], self.num[p])
            else:
                for i in range(n):
                    for p in range(self.ir[i], self.ir[i+1]):
                        yield (i, self.jc[p], 1)
        else:
            for i in range(n):
                for p in range(self.ir[i], self.ir[i+1]):
                    yield i, self.jc[p]

    def bfs(self, source):
        n = self.num_vertices()
        levels = np.ones(n, dtype=int) * -1
        parents = np.ones(n, dtype=int) * -1
        frontier = [source]
        parents[source] = source
        level = 0
        levels[source] = level

        while len(frontier) > 0:
            level += 1
            neighbors = []
            for u in frontier:
                for v in self.neighbors(u):
                    if levels[v] == -1:
                        levels[v] = level
                        parents[v] = u
                        neighbors.append(v)
            frontier = neighbors[:]

        return levels, parents

    def apsp(self):
        n = self.num_vertices()
        D = np.ones((n,n), dtype=float) * np.inf

        for i in range(n):
            D[i,i] = 0

        for i,j,w in self.edges(add_weights=True):
            D[i,j] = w

        for k in range(n):
            for i in range(n):
                for j in range(n):
                    D[i,j] = min(D[i,j], D[i,k] + D[k,j])

        return D

    def sssp(self, source):
        n = self.num_vertices()
        visited = set()
        parents = np.ones(n, dtype=int) * -1
        parents[source] = source
        costs = np.ones(n, dtype=float) * np.inf
        costs[source] = 0
        pq = []
        heap.heappush(pq, (0, source))

        while pq:
            _, u = heap.heappop(pq)
            visited.add(u)
            for v, weight in self.neighbors(u, add_weights=True):
                if v in visited:
                    continue
                newcost = costs[u] + weight
                if costs[v] > newcost:
                    parents[v] = u
                    costs[v] = newcost
                    heap.heappush(pq, (newcost, v))

        return costs, parents

def usage():
    sys.stderr.write("usage: python suite.py [graph.txt] [bfs || sssp || apsp] ([source])\n")
    sys.stderr.flush()
    sys.exit(-1)

if __name__ == "__main__":

    argc = len(sys.argv)

    if argc < 3 or argc > 4:
        usage()

    graph_filename = sys.argv[1]
    graph_algorithm = sys.argv[2]

    if argc == 4:
        source = int(sys.argv[3])

    if graph_algorithm not in ["bfs", "sssp", "apsp"]:
        usage()

    if graph_algorithm in ["bfs", "sssp"] and argc != 4:
        usage()

    if graph_algorithm == "apsp" and argc == 4:
        usage()

    G = Graph.read(graph_filename)

    if graph_algorithm == "bfs":
        levels, parents = G.bfs(source)
        sys.stdout.write(",".join(str(levels[i]) for i in range(G.num_vertices())) + '\n')
        sys.stdout.write(",".join(str(parents[i]) for i in range(G.num_vertices())) + '\n')
        sys.stdout.flush()
    elif graph_algorithm == "sssp":
        costs, parents = G.sssp(source)
        if G.weight == GraphWeight.UNWEIGHTED or G.weight == GraphWeight.WEIGHTED_INT:
            costtype = lambda x: -1 if x == np.inf else int(x)
        else:
            costtype = lambda x: -1 if x == np.inf else float(x)
        sys.stdout.write(",".join(str(costtype(costs[i])) for i in range(G.num_vertices())) + '\n')
        sys.stdout.write(",".join(str(parents[i]) for i in range(G.num_vertices())) + '\n')
        sys.stdout.flush()
    elif graph_algorithm == "apsp":
        D = G.apsp()
        for i in range(G.num_vertices()):
            sys.stdout.write(",".join(str(round(D[i,j], 5)) for j in range(G.num_vertices())) + '\n')
        sys.stdout.flush()



