import sys
import numpy as np
from enum import IntEnum

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
