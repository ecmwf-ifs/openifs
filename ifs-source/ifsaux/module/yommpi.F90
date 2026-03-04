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

MODULE YOMMPI


USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Contains identifiers used by MPI (Message Passing Interface)


INTEGER(KIND=JPIM), PARAMETER :: MINTET = 1
INTEGER(KIND=JPIM), PARAMETER :: MREALT = 2
INTEGER(KIND=JPIM), PARAMETER :: MLOGIT = 3
INTEGER(KIND=JPIM), PARAMETER :: MCHART = 4

END MODULE YOMMPI
