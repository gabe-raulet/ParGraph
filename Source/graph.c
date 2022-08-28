#include "graph.h"
#include "array.h"
#include "minheap.h"
#include "bitmap.h"
#include "infdouble.h"
#include <math.h>

void apsp(const spmat *graph, double *dist)
{
    long n, k, i, j, p;

    n = getnrows(graph);

    for (i = 0; i < n*n; ++i)
        dist[i] = -1;

    for (i = 0; i < n; ++i)
        for (p = graph->ir[i]; p < graph->ir[i+1]; ++p)
            dist[i*n + graph->jc[p]] = graph->num? graph->num[p] : 1;

    for (i = 0; i < n; ++i)
        dist[i*n + i] = 0;

    for (k = 0; k < n; ++k)
        for (i = 0; i < n; ++i)
            for (j = 0; j < n; ++j)
                dist[i*n + j] = infmin(dist[i*n + j], infplus(dist[i*n + k], dist[k*n + j]));
}

void sssp(const spmat *graph, long source, double *costs, long *parents)
{
    long n = getnrows(graph);
    bitmap *visited = bitmap_init(n);

    for (long i = 0; i < n; ++i)
    {
        parents[i] = -1;
        costs[i] = -1;
    }

    parents[source] = source;
    costs[source] = 0;

    minheap *pq = minheap_init((long)log2(n));
    minheap_insert(pq, 0, source);

    while (!minheap_empty(pq))
    {
        long u = minheap_extract(pq);
        bitmap_set(visited, u);

        for (long p = graph->ir[u]; p < graph->ir[u+1]; ++p)
        {
            long v = graph->jc[p];

            if (bitmap_get(visited, v))
                continue;

            double weight = graph->num? graph->num[p] : 1;
            double newcost = infplus(costs[u], weight);

            if (inflt(newcost, costs[v]))
            {
                parents[v] = u;
                costs[v] = newcost;
                minheap_insert(pq, newcost, v);
            }
        }
    }

    bitmap_free(visited);
    minheap_free(pq);
}

void bfs(const spmat *graph, long source, long *levels, long *parents)
{
    array_t(long) a1, a2, *frontier, *neighbors, *t;
    array_init(a1);
    array_init(a2);
    long n = getnrows(graph);

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
            for (long vi = graph->ir[u]; vi < graph->ir[u+1]; ++vi)
            {
                long v = graph->jc[vi];
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
