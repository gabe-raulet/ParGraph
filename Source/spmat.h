#pragma once

#include <stdio.h>

/* compressed sparse row */
typedef struct
{
    long m;      /* number of rows */
    long n;      /* number of columns */
    long *ir;    /* row pointers */
    long *jc;    /* column indices */
    double *num; /* value array, NULL if binary matrix */
} spmat;

#define getnrows(mat) ((mat)->m)
#define getncols(mat) ((mat)->n)
#define getnnz(mat) ((mat)->ir[getnrows(mat)])

spmat* spmat_init(long nrows, long ncols, long nnz, const long *ri, const long *ci, const double *vs);
void spmat_free(spmat *mat);
spmat* spmat_copy(const spmat *mat);
spmat* spmat_transpose(const spmat *mat);
