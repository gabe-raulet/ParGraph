#include "graph.h"
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int usage(const char *exepath)
{
    fprintf(stderr, "usage: %s [graph.txt] [bfs || sssp || apsp] ([source])\n", exepath);
    return -1;
}

void print_long_array(FILE *f, const long *arr, long n)
{
    for (long i = 0; i < n-1; ++i)
        fprintf(f, "%ld,", arr[i]);
    fprintf(f, "%ld\n", arr[n-1]);
}

void print_double_array(FILE *f, const double *arr, long n)
{
    for (long i = 0; i < n-1; ++i)
        fprintf(f, "%lg,", arr[i]);
    fprintf(f, "%lg\n", arr[n-1]);
}

int main(int argc, char *argv[])
{
    if (argc < 3 || argc > 4) return usage(*argv);

    const char *fname = argv[1];
    const char *algo = argv[2];

    int source;

    if (argc == 4) source = atoi(argv[3]);

    if (strcmp(algo, "bfs") && strcmp(algo, "sssp") && strcmp(algo, "apsp"))
        return usage(*argv);

    if ((!strcmp(algo, "bfs") || !strcmp(algo, "sssp")) && argc != 4)
        return usage(*argv);

    if (!strcmp(algo, "apsp") && argc == 4)
        return usage(*argv);

    FILE *f = fopen(fname, "r");
    spmat *g = read_graph(f);
    fclose(f);

    long n = getnrows(g);

    if (!strcmp(algo, "bfs"))
    {
        long *levels = malloc(n * sizeof(long));
        long *parents = malloc(n * sizeof(long));
        bfs(g, source, levels, parents);

        print_long_array(stdout, levels, n);
        print_long_array(stdout, parents, n);

        free(levels);
        free(parents);
    }
    else if (!strcmp(algo, "apsp"))
    {
        double *dist = malloc(n*n * sizeof(double));
        apsp(g, dist);

        for (long i = 0; i < n; ++i)
            print_double_array(stdout, &dist[i*n], n);

        free(dist);
    }

    spmat_free(g);
    return 0;
}
