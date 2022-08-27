#include "spmat.h"
#include <stdlib.h>
#include <string.h>

void spmat_free(spmat *mat)
{
    if (!mat) return;
    mat->m = mat->n = 0;
    free(mat->ir);
    free(mat->jc);
    if (mat->num) free(mat->num);
    free(mat);
}

spmat* spmat_copy(const spmat *mat)
{
    spmat *mcopy = malloc(sizeof(spmat));
    mcopy->m = mat->m;
    mcopy->n = mat->n;
    mcopy->ir = malloc((mcopy->n+1) * sizeof(long));
    mcopy->jc = malloc(getnnz(mat) * sizeof(long));
    mcopy->num = mat->num? malloc(getnnz(mat) * sizeof(double)) : NULL;
    memcpy(mcopy->ir, mat->ir, (mcopy->n+1) * sizeof(long));
    memcpy(mcopy->jc, mat->jc, getnnz(mat) * sizeof(long));
    if (mcopy->num) memcpy(mcopy->num, mat->num, getnnz(mat) * sizeof(double));
    return mcopy;
}

spmat* spmat_transpose(const spmat *mat)
{
    spmat *trmat;
    long m, n, nz, *ri, *ci, k, i, p;
    double *vs;

    m = getnrows(mat);
    n = getncols(mat);
    nz = getnnz(mat);

    ri = malloc(nz * sizeof(long));
    ci = malloc(nz * sizeof(long));
    vs = mat->num? malloc(nz * sizeof(double)) : NULL;

    for (i = k = 0; i < m; ++i)
        for (p = mat->ir[i]; p < mat->ir[i+1]; ++p)
        {
            ci[k] = i;
            ri[k] = mat->jc[p];
            if (vs) vs[k] = mat->num[p];
            ++k;
        }

    trmat = spmat_init(n, m, nz, ri, ci, vs);

    free(ri);
    free(ci);
    if (vs) free(vs);

    return trmat;
}

spmat* spmat_init(long nrows, long ncols, long nnz, const long *ri, const long *ci, const double *vs)
{
    spmat *mat;
    long *ir, *jc, *rowcnts, k, p, nzcnt;
    double *num;

    ir = calloc(nrows+1, sizeof(long));
    jc = malloc(nnz * sizeof(long));
    num = vs? malloc(nnz * sizeof(double)) : NULL;
    rowcnts = calloc(nrows, sizeof(long));

    for (k = 0; k < nnz; ++k)
        ++rowcnts[ri[k]];

    for (k = nzcnt = 0; k < nrows; ++k)
    {
        ir[k] = nzcnt;
        nzcnt += rowcnts[k];
        rowcnts[k] = ir[k];
    }

    ir[nrows] = nzcnt;

    for (k = 0; k < nnz; ++k)
    {
        p = rowcnts[ri[k]]++;
        jc[p] = ci[k];
        if (num) num[p] = vs[k];
    }

    free(rowcnts);

    mat = malloc(sizeof(spmat));
    mat->m = nrows;
    mat->n = ncols;
    mat->ir = ir;
    mat->jc = jc;
    mat->num = num;

    return mat;
}

