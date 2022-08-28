#pragma once

#include <stdint.h>

typedef struct
{
    uint64_t *v;
    long n;
} bitmap;

/**
 * @macro bitmap_items
 *
 * @param n [long] - number of bits
 *
 * @return [long] - number of longs needed to represent @n bits
 **/
#define bitmap_items(n) (((n)+63)>>6)

/**
 * @macro bitmap_set
 *
 * @param b [bitmap*] - an initialized bitmap object
 * @param i [long]    - an indexing location within @b
 *
 * Set bit @i in @b.
 **/
#define bitmap_set(b, i) ((b)->v[(i)>>6]|=(1ULL<<((i)&0x3f)))

/**
 * @macro bitmap_get
 *
 * @param b [bitmap*] - an initialized bitmap object
 * @param i [long]    - an indexing location within @b
 *
 * @return [boolean] - True if bit @i is set; False otherwise
 **/
#define bitmap_get(b, i) (((b)->v[(i)>>6]>>((i)&0x3f))&1)


/**
 * @macro bitmap_unset
 *
 * @param b [bitmap*] - an initialized bitmap object
 * @param i [long]    - an indexing location within @b
 *
 * Unset bit @i in @b.
 **/
#define bitmap_unset(b, i) ((b)->v[(i)>>6]^=(1ULL<<((i)&0x3f)))

/**
 * @macro bitmap_size
 *
 * @param b [bitmap*] - an initialized bitmap object
 *
 * @return [long] - the number of bits in @b
 **/
#define bitmap_size(b) ((b)->n)

/**
 * @func bitmap_init
 *
 * @param n [long] - number of bits
 *
 * @return [bitmap*] - an new initialized bitmap object with @n bits
 **/
bitmap* bitmap_init(long n);

/**
 * @func bitmap_free
 *
 * @param b [bitmap*] - an initialized bitmap object
 *
 * Free all the memory associated with an initialized bitmap object.
 **/
void bitmap_free(bitmap *b);

/**
 * @func bitmap_reset
 *
 * @param b [bitmap*] - an initialized bitmap object
 *
 * Unset all bits in @b.
 **/
void bitmap_reset(bitmap *b);

/**
 * @func bitmap_next
 *
 * @param b   [const bitmap*] - an initialized bitmap object
 * @param pos [long]          - position from which to find the next set bit
 *
 * @return [long] - the index of the next set bit in @b following
 **/
long bitmap_next(const bitmap *b, long pos);

/**
 * @func bitmap_indices
 *
 * @param b      [const bitmap*] - an initialized bitmap object
 * @param sparse [long*]         - a pointer to an initialized long array
 *
 * @return [long] - the size of the @sparse array
 *
 * Fills a sparse vector @sparse where each entry corresponds to
 * a set bit in @b.
 **/
long bitmap_indices(const bitmap *b, long *sparse);

/**
 * @func bitmap_numset
 *
 * @param b [const bitmap*] - an initialized bitmap object
 *
 * @return [long] - the number of bits set in @b
 **/
long bitmap_numset(const bitmap *b);

/**
 * @func bitmap_empty
 *
 * @param b [const bitmap*] - an initialized bitmap object
 *
 * @return [boolean] - True if no bits are set; False otherwise
 **/
int bitmap_empty(const bitmap *b);

/**
 * @func bitmap_union_update
 *
 * @param b [bitmap*]       - an initialized bitmap object
 * @param c [const bitmap*] - an initialized bitmap object
 *
 * Union update @b with @c.
 **/
void bitmap_union_update(bitmap *b, const bitmap *c);

/**
 * @func bitmap_apply
 *
 * @param b [const bitmap*] - an initialized bitmap object
 * @param w [long*]         - an array of longs
 * @param x [long]          - value to set
 *
 * All the indices where @b has a set bit, sets the
 * corresponding entry of @w to @x.
 **/
void bitmap_apply(const bitmap *b, long *w, long x);
