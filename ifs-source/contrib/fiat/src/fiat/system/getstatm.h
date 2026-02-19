/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


/* statm format: (/proc/self/mem under LINUX)

   in units of pages (4096 bytes)

size       total program size
resident   size of in memory portions
shared     number of the pages that are shared
trs        number of pages that are 'code'
drs        number of pages of data/stack
lrs        number of pages of library
dt         number of dirty pages

*/

struct statm
{
    int size;
    int resident;
    int shared;
    int trs;
    int drs;
    int lrs;
    int dt;
};

extern int getstatm(struct statm *sm);

