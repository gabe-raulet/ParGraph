#pragma once

typedef struct heapitem heapitem;
typedef struct minheap minheap;

minheap* minheap_init(long size);
int minheap_empty(const minheap *hp);
void minheap_insert(minheap *hp, double key, long value);
long minheap_extract(minheap *hp);
void minheap_free(minheap *hp);
