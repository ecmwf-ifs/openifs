/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>


/*
 *
 * Philippe Marguinaud, Meteo-France
 *
 */

#if defined(LINUX)

#if defined(__PGI)
#warning PGI (tested up to 19.4) does not support __sync_fetch_and_add, skipping memory_hook
/*
 * Comments/Modifications by Willem Deconinck, ECMWF
 * Problem with PGI: undefined symbol __sync_fetch_and_add
 *
 * --> Is this not GNU specific rather than LINUX?
 * + C11 has new standard API available in <stdatomic.h>,
 *   but PGI/19.4 does not yet support it.
 */
#else

static size_t align = 0;  /* Must be a multiple of sizeof (void) */
static unsigned char snan8[] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf4, 0x7f };
static unsigned char cinit[] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
static unsigned char zero0[] = { 0x00 };
static unsigned char * init = NULL;
static int sizeof_init = 0;
static int count = -1;

void * __wrap_malloc (size_t size, const void *caller)
{
  void * ptr;
  int _align = align > 0 ? align : sizeof (void*);
  int c;

  if ((c = posix_memalign (&ptr, _align, size)) != 0)
    {
      ptr = NULL; 
      printf (" c = %d, EINVAL = %d, ENOMEM = %d, align = %ld, size = %ld\n",
                c,      EINVAL,      ENOMEM,      align,       size);
    }

  if ((init != NULL) && (ptr != NULL))
    {
      unsigned char * c;
      size_t i;
      for (i = 0, c = ptr; i < size; i++)
        c[i] = init[i%sizeof_init];
    }

  if (count >= 0)
    __sync_fetch_and_add (&count, 1);

  return ptr;
}

void __attribute__((constructor))
memory_hook_init_ ()
{
  const char * MEMORY_HOOK_ALIGN = getenv ("MEMORY_HOOK_ALIGN");
  const char * MEMORY_HOOK_INIT  = getenv ("MEMORY_HOOK_INIT");
  const char * MEMORY_HOOK_COUNT = getenv ("MEMORY_HOOK_COUNT");
  if (MEMORY_HOOK_INIT)
    {
      if (strcasecmp (MEMORY_HOOK_INIT, "NAN") == 0)
        {
          init = &snan8[0];
          sizeof_init = sizeof (snan8);
        }
      else if (strcasecmp (MEMORY_HOOK_INIT, "ZERO") == 0)
        {
          init = &zero0[0];
          sizeof_init = sizeof (zero0);
        }
      else if (strncasecmp (MEMORY_HOOK_INIT, "0X", 2) == 0)
        {
          long long unsigned int x = 0;
	  const char * c;
          for (c = MEMORY_HOOK_INIT+2; *c; c++)
            if (('0' <= *c) && (*c <= '9'))
              x = 16 * x + (*c - '0');
	    else if (('a' <= *c) && (*c <= 'f'))
              x = 16 * x + (*c - 'a' + 10);
	    else if (('A' <= *c) && (*c <= 'F'))
              x = 16 * x + (*c - 'A' + 10);
            else
              break;
          init = &cinit[0];
	  memcpy (cinit, &x, sizeof (x));
          sizeof_init = sizeof (cinit);
        }
      printf (" MEMORY_HOOK_INIT = %s\n", MEMORY_HOOK_INIT);
    }
  if (MEMORY_HOOK_ALIGN)
    align = atoi (MEMORY_HOOK_ALIGN);
  if (MEMORY_HOOK_COUNT)
    count = 0;
}

void __attribute__((destructor))
memory_hook_exit_ ()
{
  if (count > 0)
    printf ("MEMORY_HOOK_COUNT = %d\n", count);
}
#endif
#endif
