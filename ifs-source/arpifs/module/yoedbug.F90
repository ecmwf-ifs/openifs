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

MODULE YOEDBUG

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     -------------------------------------------------
!*    ** *YOEDBUG* - CONTROL OPTIONS FOR DEBUGGING HELP
!     -------------------------------------------------

TYPE :: TEDBUG
INTEGER(KIND=JPIM) :: KSTPDBG(10), NSTPDBG
END TYPE TEDBUG

!!TYPE(TEDBUG), POINTER :: YREDBUG => NULL()

!     ------------------------------------------------------------------
END MODULE YOEDBUG

