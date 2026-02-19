/*
! (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
! (C) Copyright 1989- Meteo-France.
 * 
*/

#ifndef _BBT_H
#define _BBT_H


void bbt_insert_bbt(void *, void *, int (*)());
void bbt_delete_bbt(void *, void *, int (*)());

#define BBT_HEADER(self) int priority; struct self *left, *right


#endif
