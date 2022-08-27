#include "graph.h"
#include "array.h"
#include <math.h>

void bfs(const spmat *G, long source, long *levels, long *parents)
{
    array_t(long) a1, a2, *frontier, *neighbors, *t;
    array_init(a1);
    array_init(a2);
    long n = getnrows(G);

    frontier = &a1;
    neighbors = &a2;

    for (long i = 0; i < n; ++i)
        levels[i] = parents[i] = -1;

    levels[source] = 0;
    parents[source] = source;

    *array_push(*frontier) = source;

    long level = 1;

    while (!array_empty(*frontier))
    {
        array_clear(*neighbors);

        for (long ui = 0; ui < array_size(*frontier); ++ui)
        {
            long u = array_at(*frontier, ui);
            for (long vi = G->ir[u]; vi < G->ir[u+1]; ++vi)
            {
                long v = G->jc[vi];
                if (levels[v] == -1)
                {
                    levels[v] = level;
                    parents[v] = u;
                    *array_push(*neighbors) = v;
                }
            }
        }

        t = frontier;
        frontier = neighbors;
        neighbors = t;

        ++level;
    }

    array_free(a1);
    array_free(a2);
}

void write_graph(const spmat *graph, int dir, int weight, FILE *f)
{
    if (!graph || !f || weight < 0 || weight > 2 || dir < 0 || dir > 1 || getnrows(graph) != getncols(graph))
        return;

    long n = getnrows(graph);
    long m = getnnz(graph);

    if (dir == GRAPH_UNDIRECTED)
        m /= 2;

    fprintf(f, "%ld %ld %d %d\n", n, m, dir, weight);

    for (long i = 0; i < n; ++i)
        for (long p = graph->ir[i]; p < graph->ir[i+1]; ++p)
        {
            long j = graph->jc[p];

            if (dir == GRAPH_UNDIRECTED && i > j)
                continue;

            if (weight == GRAPH_UNWEIGHTED)
            {
                fprintf(f, "%ld %ld\n", i, graph->jc[p]);
            }
            else if (weight == GRAPH_WEIGHTED_INT)
            {
                long v = graph->num? ((long)ceil(graph->num[p])) : 1;
                fprintf(f, "%ld %ld %ld\n", i, j, v);
            }
            else
            {
                double v = graph->num? graph->num[p] : 1.0;
                fprintf(f, "%ld %ld %lg\n", i, j, v);
            }
        }
}

spmat* read_graph(FILE *f)
{
    spmat *graph;
    char line[1024];
    long m, n, nz, *ri, *ci, k, p, i, j;
    int asym, weight;
    double *vs, v;

    #define FREE_ALL (free(ri), free(ci), free(vs), NULL)

    if (!fgets(line, 1024, f))
        return NULL;

    if (sscanf(line, "%ld %ld %d %d", &n, &nz, &asym, &weight) != 4)
        return NULL;

    m = nz;

    if (!asym) m *= 2;

    ri = malloc(m * sizeof(long));
    ci = malloc(m * sizeof(long));

    vs = (weight != GRAPH_UNWEIGHTED)? malloc(m * sizeof(double)) : NULL;

    for (k = p = 0; k < nz; ++k)
    {
        if (!fgets(line, 1024, f)) return FREE_ALL;
        if (!vs && (sscanf(line, "%ld %ld", &i, &j) != 2)) return FREE_ALL;
        else if (vs && (sscanf(line, "%ld %ld %lg", &i, &j, &v) != 3)) return FREE_ALL;

        ri[p] = i;
        ci[p] = j;
        if (vs) vs[p] = v;

        ++p;

        if (!asym)
        {
            ri[p] = j;
            ci[p] = i;
            if (vs) vs[p] = v;
            ++p;
        }
    }

    graph = spmat_init(n, n, m, ri, ci, vs);

    FREE_ALL;

    return graph;
}
