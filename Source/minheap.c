#include "minheap.h"
#include "array.h"
#include "infdouble.h"
#include <math.h>
#include <limits.h>

struct heapitem
{
    double key;
    long value;
};

struct minheap
{
    array_t(heapitem) heap;
};

static inline unsigned long log2ul(unsigned long val)
{
    if (!val) return ULONG_MAX;
    if (val == 1) return 0;
    unsigned long r = 0;
    while (val > 1) { val >>= 1; ++r; }
    return r;
}

#define MIN(a, b) (((a) < (b))? (a) : (b))
#define MAX(a, b) (((a) > (b))? (a) : (b))

#define LEFT_CHILD(i) (2*(i) + 1)
#define RIGHT_CHILD(i) (2*(i) + 2)
#define PARENT(i) (((i)-1)/2)

minheap* minheap_init(long size)
{
    minheap *hp = malloc(sizeof(minheap));
    array_init(hp->heap);
    array_reserve(hp->heap, size);
    return hp;
}

void minheap_free(minheap *hp)
{
    if (!hp) return;
    array_free(hp->heap);
    free(hp);
}

int minheap_empty(const minheap *hp)
{
    return array_empty(hp->heap);
}

void minheap_insert(minheap *hp, double key, long value)
{
    long i = array_size(hp->heap);
    *array_push(hp->heap) = (heapitem){key, value};

    heapitem *heap = &array_at(hp->heap, 0);
    heapitem t;

    while (i > 0 && inflt(heap[i].key, heap[PARENT(i)].key))
    {
        t = heap[PARENT(i)];
        heap[PARENT(i)] = heap[i];
        heap[i] = t;
        i = PARENT(i);
    }
}

long minheap_extract(minheap *hp)
{
    heapitem *heap = &array_at(hp->heap, 0);
    heapitem t;
    long minitem = heap[0].value;

    heap[0] = array_pop(hp->heap);
    long size = array_size(hp->heap);

    long i = 0;
    double mykey, leftkey, rightkey;
    int left_lt_me, right_lt_me, left_lt_right;

    while (LEFT_CHILD(i) < size)
    {
        mykey = heap[i].key;
        leftkey = heap[LEFT_CHILD(i)].key;

        left_lt_me = inflt(leftkey, mykey);

        if (RIGHT_CHILD(i) < size)
        {
            rightkey = heap[RIGHT_CHILD(i)].key;
            right_lt_me = inflt(rightkey, mykey);
            left_lt_right = inflt(leftkey, rightkey);
        }
        else left_lt_right = 1;

        if (left_lt_me || right_lt_me)
        {
            if (left_lt_right)
            {
                t = heap[i];
                heap[i] = heap[LEFT_CHILD(i)];
                heap[LEFT_CHILD(i)] = t;
                i = LEFT_CHILD(i);
            }
            else
            {
                t = heap[i];
                heap[i] = heap[RIGHT_CHILD(i)];
                heap[RIGHT_CHILD(i)] = t;
                i = RIGHT_CHILD(i);
            }
        }
        else break;
    }
    return minitem;
}
