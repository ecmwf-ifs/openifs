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

MODULE YOMPGOM

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMPGO,   ONLY : NLENCHA

IMPLICIT NONE

SAVE

!*  DECK FOR CODING DDH OUTPUT IN PSEUDO-GRIB

! CHAPG  - CHARACTER BLOCK
! REAPG(NLENREA)  - REAL BLOCK
! NINDPG(NLENIND) - INDEX ARRAY
! NINTPG(NLENINT) - INTEGER BLOCK

CHARACTER(LEN=NLENCHA) :: CHAPG

REAL(KIND=JPRB),ALLOCATABLE:: REAPG(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NINDPG(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NINTPG(:)

END MODULE YOMPGOM
