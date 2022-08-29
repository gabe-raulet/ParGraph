#pragma once

#include "spmat.h"

#define GRAPH_UNDIRECTED (0)
#define GRAPH_DIRECTED   (1)
#define GRAPH_UNWEIGHTED (0)
#define GRAPH_WEIGHTED_INT (1)
#define GRAPH_WEIGHTED_REAL (2)

spmat* read_graph(FILE *f);
void write_graph(const spmat *graph, int dir, int weight, FILE *f);

void bfs(const spmat *graph, long source, long *levels, long *parents, long limit);
void uy(const spmat *graph, long source, double constant, long *levels);
void sssp(const spmat *graph, long source, double *costs, long *parents);
void apsp(const spmat *graph, double *dist);
