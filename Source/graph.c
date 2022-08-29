#include "graph.h"
#include "array.h"
#include "minheap.h"
#include "bitmap.h"
#include "infdouble.h"
#include <assert.h>
#include <math.h>
#include <omp.h>

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
        costs[i] = -1;
        if (parents) parents[i] = -1;
    }

    costs[source] = 0;
    if (parents) parents[source] = source;

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
                costs[v] = newcost;
                if (parents) parents[v] = u;
                minheap_insert(pq, newcost, v);
            }
        }
    }

    bitmap_free(visited);
    minheap_free(pq);
}

/**
 * @func fisher_yates
 *
 * @param arr [long*] - array of >= k longs where random sample goes
 * @param k   [long]  - number of integers to sample
 * @param n   [long]  - sample from [0..n-1]
 *
 * Randomly sample k integers without replacement from [0..n-1].
 **/
void fisher_yates(long *arr, long k, long n)
{
    long *x, l, i, r, t;

    x = malloc(n * sizeof(long));

    for (l = 0; l < n; ++l)
        x[l] = l;

    for (i = 0; i < k; ++i)
    {
        r = rand() % (n-i);
        t = x[r];
        x[r] = x[n-i-1];
        x[n-i-1] = t;
    }

    for (i = 0; i < k; ++i)
        arr[i] = x[n-k+i];

    free(x);
}

void uy_seq(const spmat *graph, long source, double constant, long *levels)
{
    long n;  /* number of vertices */
    long nS; /* number of distinguished vertices */
    long *S; /* array of distinguished vertices */

    n = getnrows(graph);
    nS = (long)floor(constant*sqrt(n)*log2(n));
    S = malloc((nS+1) * sizeof(long));

    fisher_yates(S, nS, n);

    long t = -1;

    for (long i = 0; i < nS; ++i)
        if (S[i] == source)
        {
            t = S[0];
            S[0] = S[i];
            S[i] = t;
            break;
        }

    if (t < 0)
    {
        S[nS++] = S[0];
        S[0] = source;
    }

    long **S_levels = malloc(nS * sizeof(long*));

    for (long i = 0; i < nS; ++i)
        S_levels[i] = malloc(n * sizeof(long));

    long limit = (long)ceil(sqrt(n));

    for (long i = 0; i < nS; ++i)
    {
        long u = S[i];
        bfs(graph, u, S_levels[i], NULL, limit);
    }

    array_t(long) ri_H, ci_H;
    array_t(double) vs_H;

    array_init(ri_H);
    array_init(ci_H);
    array_init(vs_H);

    for (long xH = 0; xH < nS; ++xH)
        for (long yH = 0; yH < nS; ++yH)
        {
            if (xH == yH) continue;

            long yG = S[yH];
            double xycost = S_levels[xH][yG];

            if (xycost >= 0)
            {
                *array_push(ri_H) = xH;
                *array_push(ci_H) = yH;
                *array_push(vs_H) = xycost;
            }
        }

    long mH = array_size(ri_H);
    long *ri = array_release(ri_H);
    long *ci = array_release(ci_H);
    double *vs = array_release(vs_H);

    spmat *H = spmat_init(nS, nS, mH, ri, ci, vs);

    free(ri);
    free(ci);
    free(vs);

    double **costs = malloc(nS * sizeof(double*));

    for (long i = 0; i < nS; ++i)
        costs[i] = malloc(nS * sizeof(double));

    for (long s = 0; s < nS; ++s)
        sssp(H, s, costs[s], NULL);

    spmat_free(H);

    for (long i = 0; i < n; ++i)
        levels[i] = -1;

    levels[source] = 0;

    for (long v = 0; v < n; ++v)
    {
        double minpath = -1;

        for (long x = 0; x < nS; ++x)
        {
            double Psx = costs[0][x];
            double Pxv = (double)S_levels[x][v];
            double Psv = infplus(Psx, Pxv);

            minpath = infmin(minpath, Psv);
        }

        levels[v] = (long)minpath;
    }

    free(S);

    for (long i = 0; i < nS; ++i)
    {
        free(S_levels[i]);
        free(costs[i]);
    }

    free(S_levels);
    free(costs);
}

void uy(const spmat *graph, long source, double constant, long *levels)
{
    long n;  /* number of vertices */
    long nS; /* number of distinguished vertices */
    long *S; /* array of distinguished vertices */

    n = getnrows(graph);
    nS = (long)floor(constant*sqrt(n)*log2(n));
    S = malloc((nS+1) * sizeof(long));

    fisher_yates(S, nS, n);

    long t = -1;

    for (long i = 0; i < nS; ++i)
        if (S[i] == source)
        {
            t = S[0];
            S[0] = S[i];
            S[i] = t;
            break;
        }

    if (t < 0)
    {
        S[nS++] = S[0];
        S[0] = source;
    }

    int nthreads = 1;

    #pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }

    long **S_levels = malloc(nS * sizeof(long*));

    for (long i = 0; i < nS; ++i)
        S_levels[i] = malloc(n * sizeof(long));

    long limit = (long)ceil(sqrt(n));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthrds = omp_get_num_threads();

        assert(nthrds == nthreads);

        long chunksize = nS / nthreads;
        long start = chunksize * ((long)tid);
        long end = (tid != nthreads-1)? start + chunksize : nS;

        for (long i = start; i < end; ++i)
        {
            long u = S[i];
            bfs(graph, u, S_levels[i], NULL, limit);
        }
    }

    array_t(long) ri_H, ci_H;
    array_t(double) vs_H;

    array_init(ri_H);
    array_init(ci_H);
    array_init(vs_H);

    for (long xH = 0; xH < nS; ++xH)
        for (long yH = 0; yH < nS; ++yH)
        {
            if (xH == yH) continue;

            long yG = S[yH];
            double xycost = S_levels[xH][yG];

            if (xycost >= 0)
            {
                *array_push(ri_H) = xH;
                *array_push(ci_H) = yH;
                *array_push(vs_H) = xycost;
            }
        }

    long mH = array_size(ri_H);
    long *ri = array_release(ri_H);
    long *ci = array_release(ci_H);
    double *vs = array_release(vs_H);

    spmat *H = spmat_init(nS, nS, mH, ri, ci, vs);

    free(ri);
    free(ci);
    free(vs);

    double **costs = malloc(nS * sizeof(double*));

    for (long i = 0; i < nS; ++i)
        costs[i] = malloc(nS * sizeof(double));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthrds = omp_get_num_threads();

        assert(nthrds == nthreads);

        long chunksize = nS / nthreads;
        long start = chunksize * ((long)tid);
        long end = (tid != nthreads-1)? start + chunksize : nS;

        for (long s = start; s < end; ++s)
            sssp(H, s, costs[s], NULL);
    }

    spmat_free(H);

    #pragma omp parallel for schedule(static, 1)
    for (long i = 0; i < n; ++i)
        levels[i] = -1;

    levels[source] = 0;

    #pragma omp parallel for schedule(static, 1)
    for (long v = 0; v < n; ++v)
    {
        double minpath = -1;

        for (long x = 0; x < nS; ++x)
        {
            double Psx = costs[0][x];
            double Pxv = (double)S_levels[x][v];
            double Psv = infplus(Psx, Pxv);

            minpath = infmin(minpath, Psv);
        }

        levels[v] = (long)minpath;
    }

    free(S);

    for (long i = 0; i < nS; ++i)
    {
        free(S_levels[i]);
        free(costs[i]);
    }

    free(S_levels);
    free(costs);
}

void bfs(const spmat *graph, long source, long *levels, long *parents, long limit)
{
    array_t(long) a1, a2, *frontier, *neighbors, *t;
    array_init(a1);
    array_init(a2);
    long n = getnrows(graph);

    frontier = &a1;
    neighbors = &a2;

    for (long i = 0; i < n; ++i)
    {
        levels[i] = -1;
        if (parents) parents[i] = -1;
    }

    levels[source] = 0;
    if (parents) parents[source] = source;

    *array_push(*frontier) = source;

    long level = 0;

    if (limit <= 0) limit = n+1;

    while (!array_empty(*frontier))
    {

        if (level++ >= limit)
            break;

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
                    if (parents) parents[v] = u;
                    *array_push(*neighbors) = v;
                }
            }
        }

        t = frontier;
        frontier = neighbors;
        neighbors = t;
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
