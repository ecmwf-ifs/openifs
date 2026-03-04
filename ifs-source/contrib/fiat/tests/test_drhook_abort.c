/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/*
 * This program initializes DR_HOOK and aborts in a DR_HOOK-traced callstack with a call to
 * `drhook_abort(...)` which delegates to user-customizable abort function.
 * By default the abort function is ABOR1. This test can be configured to instead set a
 * custom abort function by setting in the environment "CUSTOM_ABORT=1"
 * Unless "DR_HOOK_IGNORE_SIGNALS=6" is specified in the environment, DR_HOOK will also catch
 * the SIGABRT and print the drhook_calltree when "CUSTOM_ABORT=1"
 * 
 * Note that the abort may be disabled altogether with environment "ABORT=0", which is useful
 * to see if the program manages to succeed otherwise.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "drhook.h"
#include "mpl.h"

int do_abort();
int mpl();

void custom_abort( const char* file, int line, const char* txt ) {
    fprintf(stderr,"custom_abort(%s,%d,\"%s\")\n",file,line,txt);
    abort();
}

void function_2 () {
    DRHOOK_START(function_2);
    static int count = 0;
    ++count;
    usleep(100); 

    if( count == 3 ) {
        if( do_abort() ) {
            drhook_abort(__FILE__,__LINE__,"abort from call2");
        }
    }
    DRHOOK_END();
}

void function_1 () {
    int i;
    DRHOOK_START(function_1);
    for( i=0; i<5; ++i ) {
        function_2();
    }
    DRHOOK_END();
}

void setup_test(int argc, char* argv[]) {
    //ec_args(argc,argv);

    const int OVERWRITE      = 1;
    const int DONT_OVERWRITE = 0;
    setenv("DR_HOOK",        "1", OVERWRITE );
    setenv("DR_HOOK_SILENT", "1", DONT_OVERWRITE );

    char* env = getenv("CUSTOM_ABORT");
    int use_custom_abort = env ? atoi(env) : 0;
    if( use_custom_abort ) {
        drhook_set_abort( custom_abort );
    }
    if( mpl() ) {
        mpl_init();
    }

    drhook_init(argc,argv);
}
void end_test() {
    if( mpl() ) {
        mpl_end();
    }
}

int do_abort() {
    char* env = getenv("ABORT");
    return env ? atoi(env) : 1;
}

int mpl() {
    char* env = getenv("MPL");
    return env ? atoi(env) : 0;
}

int main(int argc, char* argv[]) {
    setup_test(argc,argv);
    DRHOOK_START(main);
    function_1();
    DRHOOK_END();
    end_test();
    printf("test completed\n");
}
