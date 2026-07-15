/*
 * (C) Copyright 2022- ECMWF.
 * (C) Copyright 2022- Meteo-France.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUDDY_ALLOC_IMPLEMENTATION
#include "buddy_alloc.h"
#undef BUDDY_ALLOC_IMPLEMENTATION

typedef struct
{
  void * metadata;
  void * arena;
  struct buddy * buddy;
} buddy_t;

static buddy_t B;

static void dev_alloc_init ()
{
  static int done = 0;
  if (done)
    return;
  done++;


  long long int dev_alloc_size = 64LL * 1024LL * 1024LL * 1024LL; // 64Gb
  char * DEV_ALLOC_SIZE = getenv ("DEV_ALLOC_SIZE");
  if (DEV_ALLOC_SIZE)
    sscanf (DEV_ALLOC_SIZE, "%lld", &dev_alloc_size);

  printf (" DEV_ALLOC_SIZE = %lld\n", dev_alloc_size);

  B.metadata = malloc (buddy_sizeof (dev_alloc_size));
  B.arena    = malloc (dev_alloc_size);
  B.buddy    = buddy_init (B.metadata, B.arena, dev_alloc_size);

}

static void __attribute__((destructor)) dev_alloc_exit ()
{
  if (B.metadata)
    free (B.metadata); 
  B.metadata = NULL;
  if (B.arena)
    free (B.arena);    
  B.arena = NULL;
  B.buddy = NULL;
}

void dev_malloc (size_t siz, void ** ptr)
{
  dev_alloc_init ();
  *ptr = buddy_malloc (B.buddy, siz);
  if (*ptr == NULL)
    abort ();
}

void dev_free (void * ptr)
{
  dev_alloc_init ();
  buddy_free (B.buddy, ptr);
}

