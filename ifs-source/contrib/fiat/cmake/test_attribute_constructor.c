/*
 * (C) Copyright 2022- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>

static int constructor_called = 0;

void __attribute__((constructor)) constructor() {
  printf("constructor()\n");
  constructor_called = 1;
}

int main() {
  printf("main()\n");
  if (constructor_called) {
    return 0; // success
  }
  return 1; // error
}

