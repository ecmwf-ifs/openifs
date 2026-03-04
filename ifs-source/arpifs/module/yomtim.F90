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

MODULE YOMTIM

USE PARKIND1  ,ONLY : JPRD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Time and date of the real world (both official and computer)

! RSTART : CPU time between the start of the job and last call to the timer
! RVSTART : vector CPU time between the start of the job and last call to the
! RTIMEF : wallclock time between the last two calls to the timer

REAL(KIND=JPRD) :: RSTART
REAL(KIND=JPRD) :: RVSTART
REAL(KIND=JPRD) :: RTIMEF

!     ------------------------------------------------------------------
END MODULE YOMTIM
