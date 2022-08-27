#include "graph.h"
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[])
{
    int dir = atoi(argv[2]);
    int weight = atoi(argv[3]);

    FILE *f = fopen(argv[1], "r");
    spmat *g = read_graph(f);
    fclose(f);

    write_graph(g, dir, weight, stdout);

    return 0;
}
