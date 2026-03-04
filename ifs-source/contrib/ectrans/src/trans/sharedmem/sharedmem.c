/*
 * (C) Copyright 2015- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <stdlib.h>

void sharedmem_malloc_bytes (void** ptr, size_t bytes)
{ 
  *ptr = malloc(bytes);
}

void sharedmem_free(void** ptr)
{ 
  free(*ptr);
}

void sharedmem_advance_bytes  (void** ptr, size_t bytes)
{ 
  char** char_ptr = (char**)ptr;
  *char_ptr += bytes;
}
