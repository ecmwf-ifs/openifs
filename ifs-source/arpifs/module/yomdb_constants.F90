! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMDB_CONSTANTS

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

PRIVATE :: JPIM

INTEGER(KIND=JPIM), PARAMETER :: MDIDB = 2147483647
INTEGER(KIND=JPIM), PARAMETER :: JPMAXSIMVIEWS = 24
INTEGER(KIND=JPIM), PARAMETER :: JPMAX_AUX_CASES = 5 ! How many ROBAUX<n> : ROBAUX1, ROBAUX2, ..., ROBAUX5
INTEGER(KIND=JPIM), PARAMETER :: JP_MAXDESC = 10
INTEGER(KIND=JPIM), PARAMETER :: JP_MAXDUPL = 16
INTEGER(KIND=JPIM), PARAMETER :: JP_NUMAUX   = 9 ! Must be the same as $NUMAUX in odb/ddl/cma.h
INTEGER(KIND=JPIM), PARAMETER :: JP_NUMDIAG  = 1 ! Must be the same as $NUMDIAG in odb/ddl/cma.h
INTEGER(KIND=JPIM), PARAMETER :: JP_NUMTHBOX = 3 ! Must be the same as $NUMTHBOX in odb/ddl/cma.h

!***  Database handles
INTEGER(KIND=JPIM),PARAMETER :: JPMAXSIMDB = 5   ! Max. no. of simultaneously opened DBs

END MODULE YOMDB_CONSTANTS
