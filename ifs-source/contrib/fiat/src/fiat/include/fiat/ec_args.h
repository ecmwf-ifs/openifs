/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/* ec_args.h 
 *
 * Contains routines to store and retrieve command-line arguments as presented via C-main: "int main(int argc, char* argv[])"
 * A Fortran interface is available in ec_args_mod.F90
 * The `ec_args(...)` routine should be called as early as possible
 * 
 * Author:   Sami Saarinen, ECMWF, 27-Apr-2006
 * Modified: Willem Deconinck, ECMWF, 1-Jul-2021
 * 
 * An example C program:
 *
 *      #include "ec_args.h"  
 *      int main( int argc, char* argv[] ) {
 *          ec_args(argc,argv);
 *          int num_args     = ec_argc();
 *          const char* name = ec_argv()[0];
 *      }
 * 
 */

#ifndef _EC_ARGS_H_
#define _EC_ARGS_H_

#if defined(__cplusplus)
extern "C" {
#endif

/* Register command-line arguments as presented by C-main : "int main(int argc, char** argv)" */
void ec_args(int argc, char* argv[]);

/* Return number of arguments, including the program-name, as presented by C-main */
int ec_argc(void);

/* Return arguments (0-terminated), as presented by C-main */ 
const char* const* ec_argv(void);

#if defined(__cplusplus)
} // extern "C"
#endif

#endif /* _EC_ARGS_H_ */

