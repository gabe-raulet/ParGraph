#pragma once

#include "spmat.h"

#define GRAPH_UNDIRECTED (0)
#define GRAPH_DIRECTED   (1)
#define GRAPH_UNWEIGHTED (0)
#define GRAPH_WEIGHTED_INT (1)
#define GRAPH_WEIGHTED_REAL (2)

spmat* read_graph(FILE *f);
void write_graph(const spmat *graph, int dir, int weight, FILE *f);

void bfs(const spmat *graph, long source, long *levels, long *parents);
void sssp(const spmat *graph, long source, long *costs, long *parents);
void apsp(const spmat *graph, double *dist);