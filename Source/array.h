#pragma once

#include <stdlib.h>

/**
 * @type array_t(type)
 *
 * A resizeable generic vector for templated types.
 **/
#define array_t(type) struct { type *a; long l, m; }

/**
 * @macro array_init
 *
 * @param v [array_t] - an unintialized array_t object
 *
 * Initialize a new array_t object.
 **/
#define array_init(v) ((v).l=(v).m=0, (v).a=0)

/**
 * @macro array_free
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Free all the memory associated with @v.
 **/
#define array_free(v) (free((v).a))

/**
 * @macro array_at
 *
 * @param v [array_t] - an initialized array_t object
 * @param i [long]    - an index into @v
 *
 * Retrieves the @ith item of @v. Undefined if @i is an
 * invalid index.
 **/
#define array_at(v, i) ((v).a[(i)])

/**
 * @macro array_size
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Return the number of elements currently stored in @v.
 **/
#define array_size(v) ((v).l)

/**
 * @macro array_empty
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Return True if @v has no elements and False otherwise.
 **/
#define array_empty(v) (!(v).l)

/**
 * @macro array_pop
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Remove the top element of @v and then return it.
 * Undefined if @v is empty.
 **/
#define array_pop(v) ((v).a[--(v).l])

/**
 * @macro array_clear
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Remove all the elements from @v.
 **/
#define array_clear(v) ((v).l=0)

/**
 * @macro array_reserve
 *
 * @param v [array_t] - an initialized array_t object
 * @param n [long]    - a scalar
 *
 * Allocate enough memory to @v so that at least @n items
 * can be added without another memory allocation required.
 **/
#define array_reserve(v, n) (((v).m=( ((n)>=1)? (n) : 1)), (v).a=realloc((v).a, sizeof(*(v).a)*(v).m))

/**
 * @macro array_push
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Push an item onto the end of @v. The usage is a bit unusual.
 * In order to add the item @item onto the end of @v, we would
 * do the following command:
 *
 *     *array_push(v) = item;
 **/
#define array_push(v) ((((v).l==(v).m)? (((v).m=((v).m? 2*(v).m : 2)), (v).a=realloc((v).a, sizeof(*(v).a)*(v).m)) : 0), ((v).a + ((v).l++)))

/**
 * @macro array_release
 *
 * @param v [array_t] - an initialized array_t object
 *
 * Release the memory owned by @v and return as a new standalone array.
 **/
#define array_release(v) ((((v).m>(v).l)? ((v).a=realloc((v).a, sizeof(*(v).a)*(v).l)) : 0), (v).l=(v).m=0, (v).a)
