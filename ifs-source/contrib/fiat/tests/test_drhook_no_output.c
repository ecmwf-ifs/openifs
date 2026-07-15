/*
 * (C) Copyright 2003- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdlib.h>
#include "drhook.h"

const int OVERWRITE      = 1;
const int DONT_OVERWRITE = 0;

void function_2 () {
    DRHOOK_START();
    DRHOOK_END();
}

void function_1 () {
    int i;
    DRHOOK_START();
    for( i=0; i<2; ++i ) {
        function_2();
    }
    DRHOOK_END();
}

void setup_test(int argc, char* argv[]) {
    setenv("DR_HOOK",       "1", OVERWRITE);
    setenv("DR_HOOK_SILENT","1", DONT_OVERWRITE);
    drhook_init(argc,argv);

    if( ! drhook_active() ) {
        drhook_abort(__FILE__,__LINE__,"drhook is supposed to be activated");
    }
}

int main(int argc, char* argv[]) {
    setup_test(argc,argv);

    DRHOOK_START();
    function_1();
    DRHOOK_END();
}
