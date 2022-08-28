#include "bitmap.h"
#include "array.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

bitmap* bitmap_init(long n)
{
    long nitems;
    bitmap *b;

    assert(sizeof(uint64_t) == 8);

    b = malloc(sizeof(bitmap));
    nitems = bitmap_items(n);
    b->v = calloc(nitems, 8);
    b->n = n;

    return b;
}

void bitmap_free(bitmap *b)
{
    b->n = 0;
    free(b->v);
    free(b);
}

void bitmap_reset(bitmap *b)
{
    memset(b->v, 0, bitmap_items(b->n) * sizeof(uint64_t));
}

long bitmap_next(const bitmap *b, long pos)
{
    long next_bit;
    int bit_offset;
    uint64_t *it, t, *end;

    next_bit = pos;
    bit_offset = pos&0x3f;
    it = b->v + (pos>>6);
    t = (bit_offset != 63)? *it >> (bit_offset+1) : 0;
    end = b->v + bitmap_items(b->n);

    if (!t)
    {
        next_bit &= 0xffffffc0;
        while (!t)
        {
            if (++it >= end) return -1;
            t = *it;
            next_bit += 64;
        }
    }
    else { ++next_bit; }

    while (!(t&1)) { t >>= 1; ++next_bit; }

    return next_bit;
}

long bitmap_indices(const bitmap *b, long *sparse)
{
    long pos, i;

    i = pos = 0;

    if (bitmap_get(b, pos))
        sparse[i++] = pos;

    while (1)
    {
        if ((pos = bitmap_next(b, pos)) != -1)
            sparse[i++] = pos;
        else break;
    }

    return i;
}

long bitmap_numset(const bitmap *b)
{
    uint64_t n;
    long nitems, numset;

    nitems = bitmap_items(b->n);
    numset = 0;

    for (long i = 0; i < nitems; ++i)
    {
        n = b->v[i];

        while (n)
        {
            n = n & (n - 1);
            ++numset;
        }
    }

    return numset;
}

int bitmap_empty(const bitmap *b)
{
    long nitems = bitmap_items(b->n);

    for (long i = 0; i < nitems; ++i)
        if (b->v[i] != 0)
            return 0;

    return 1;
}

void bitmap_union_update(bitmap *b, const bitmap *c)
{
    long nitems, i;

    if (bitmap_size(b) != bitmap_size(c))
        return;

    nitems = bitmap_items(b->n);

    for (i = 0; i < nitems; ++i)
        b->v[i] |= c->v[i];
}

void bitmap_apply(const bitmap *b, long *w, long x)
{
    long pos = 0;

    if (bitmap_get(b, pos))
        w[pos] = x;

    while (1)
    {
        if ((pos = bitmap_next(b, pos)) != -1)
            w[pos] = x;
        else break;
    }
}

void bitmap_print(const bitmap *b)
{
    long n = b->n;

    char *s = malloc(n+1);
    s[n] = 0;

    for (long i = 0; i < n; ++i)
        s[i] = '0';

    long pos = 0;

    if (bitmap_get(b, pos))
        s[pos] = '1';

    while (1)
    {
        if ((pos = bitmap_next(b, pos)) != -1)
            s[pos] = '1';
        else break;
    }

    printf("%s\n", s);
    free(s);
}
