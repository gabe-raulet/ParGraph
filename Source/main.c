#include "graph.h"
#include "spmat.h"
#include "minheap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static int usage(const char *exepath)
{
    fprintf(stderr, "usage: %s [graph.txt] [bfs || sssp || apsp || uy] ([source]) ([constant])\n", exepath);
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
    srand(time(0)*getpid());

    if (argc < 3 || argc > 5) return usage(*argv);

    const char *fname = argv[1];
    const char *algo = argv[2];

    int source;
    double constant;

    if (argc == 4 || argc == 5) source = atoi(argv[3]);

    if (argc == 5) constant = atof(argv[4]);

    if (strcmp(algo, "bfs") && strcmp(algo, "sssp") && strcmp(algo, "apsp") && strcmp(algo, "uy"))
        return usage(*argv);

    if ((!strcmp(algo, "bfs") || !strcmp(algo, "sssp")) && argc != 4)
        return usage(*argv);

    if (!strcmp(algo, "uy") && argc != 5)
        return usage(*argv);

    if (!strcmp(algo, "apsp") && argc >= 4)
        return usage(*argv);

    FILE *f = fopen(fname, "r");
    spmat *g = read_graph(f);
    fclose(f);

    long n = getnrows(g);

    if (!strcmp(algo, "bfs"))
    {
        long *levels = malloc(n * sizeof(long));
        long *parents = malloc(n * sizeof(long));
        bfs(g, source, levels, parents, 0);

        print_long_array(stdout, levels, n);
        print_long_array(stdout, parents, n);

        free(levels);
        free(parents);
    }
    else if (!strcmp(algo, "uy"))
    {
        long *levels = malloc(n * sizeof(long));
        uy(g, source, constant, levels);

        print_long_array(stdout, levels, n);

        free(levels);
    }
    else if (!strcmp(algo, "apsp"))
    {
        double *dist = malloc(n*n * sizeof(double));
        apsp(g, dist);

        for (long i = 0; i < n; ++i)
            print_double_array(stdout, &dist[i*n], n);

        free(dist);
    }
    else /* sssp */
    {
        double *costs = malloc(n * sizeof(double));
        long *parents = malloc(n * sizeof(long));

        sssp(g, source, costs, parents);

        print_double_array(stdout, costs, n);
        print_long_array(stdout, parents, n);


    }

    spmat_free(g);
    return 0;
}
