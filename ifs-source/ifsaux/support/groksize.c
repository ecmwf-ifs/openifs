/*
 * (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
 * (C) Copyright 1989- Meteo-France.
 * 
*/


/*
 * This routine is used to find out the size(in bytes) of a fortran entity.
 * Example : 
 * TYPE(TOTO) :: T(2)
 * CALL GROKSIZE (KSIZE, T(1), T(2))
 */
void groksize_ (int * ksize, const char * c1, const char * c2)
{
  *ksize = c1 - c2;
}
