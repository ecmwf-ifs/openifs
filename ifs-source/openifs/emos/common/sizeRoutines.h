/**
* Copyright 1981- ECMWF.
*
* This software is licensed under the terms of the Apache Licence 
* Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*
* In applying this licence, ECMWF does not waive the privileges and immunities 
* granted to it by virtue of its status as an intergovernmental organisation 
* nor does it submit to any jurisdiction.
*/
#ifndef SIZE_ROUTINES_H
#define SIZE_ROUTINES_H

static fortint prodsize(fortint, char *, fortint, fortint *,
                        fortint (*fileRead)(), void *);
static fortint gribsize(char * , fortint, fortint * , 
                        fortint (*fileRead)(), void * );
static fortint bufrsize(char * , fortint, fortint * , 
                        fortint (*fileRead)(), void * );
static fortint tide_budg_size(char *, fortint, fortint *,
                        fortint (*fileRead)(), void *);
static fortint lentotal(char *, fortint *, fortint, fortint , fortint ,
                        fortint , fortint (*fileRead)(), void *);
static fortint waveLength(char *);

static fortint crex_size( void * );

#endif /* end of  SIZE_ROUTINES_H */
