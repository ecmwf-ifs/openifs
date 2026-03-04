/*
 * (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
 * (C) Copyright 1989- Meteo-France.
 * 
*/

#include <stdlib.h>

#include "bbt.h"

typedef struct bbt_treap {
    BBT_HEADER(bbt_treap);
} bbt;

static int pseudo_random(void) {
static int x0=5341;

    x0 = (22611*x0 + 10) % 44071;
    return x0;
}

static bbt *rotate_right(bbt *t) {
bbt *temp;

    temp = t->left;
    t->left = t->left->right;
    temp->right = t;

    return temp;
}

static bbt *rotate_left(bbt *t) {
bbt *temp;

    temp =  t->right;
    t->right = t->right->left;
    temp->left = t;

    return temp;
}

static bbt *delete_root(bbt *t) {
bbt *temp;

    if (t->left == NULL)
	return t->right;

    if (t->right == NULL)
	return t->left;

    if (t->left->priority > t->right->priority) {
	temp = rotate_right(t);
	temp->right = delete_root(t);

    } else {
	temp = rotate_left(t);
	temp->left = delete_root(t);
    }

    return temp;
}

static bbt *delete_treap(bbt *old, bbt *t, int (*compare)()) {
int c;

    if (t == NULL)
	return NULL; 

    c = compare(old, t);

    if (c == 0)
	t = delete_root(t);

    else if (c < 0)
	t->left = delete_treap(old, t->left, compare);

    else if (c > 0)
	t->right = delete_treap(old, t->right, compare);

    return t;
}

void bbt_delete_bbt(void *root, void *old, int (*compare)()) {
bbt **t;

    t = (bbt **) root; 

    *t = delete_treap((bbt *) old, *t, compare);
}

static bbt *insert(bbt *new, bbt *t, int (*compare)()) {
int c;

    if (t == NULL)
	return new;

    c = (*compare)(new, t);

    if (c == 0)
	return NULL;

    if (c < 0) {
	t->left = insert(new, t->left, compare);
	if (t->priority < t->left->priority)
	    t = rotate_right(t);

    } else {
	t->right = insert(new, t->right, compare);
	if (t->priority < t->right->priority)
	    t = rotate_left(t);
    }


  return t;
}

void bbt_insert_bbt(void *root, void *new, int (*compare)()) {
bbt **r, *n;

    r = (bbt **) root; 
    n = (bbt *) new;

    n->priority = pseudo_random();
    *r = insert(n, *r, compare);
}

