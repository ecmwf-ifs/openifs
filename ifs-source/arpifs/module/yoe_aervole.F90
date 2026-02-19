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

MODULE YOE_AERVOLE

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOE_AERVOLE* - CONTROL PARAMETERS FOR ERUPTING VOLCANO
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: RVEMHGHT(19)
REAL(KIND=JPRB) :: RVEMMAS(19,160)
REAL(KIND=JPRB) :: RVMASSEM(19)

REAL(KIND=JPRB) :: RVMASSVI(0:19)
REAL(KIND=JPRB) :: RVHGHTEM(0:19)

INTEGER(KIND=JPIM) :: JVEMTIM(160)
INTEGER(KIND=JPIM) :: NVDATES
!------------------------------------------------------------------------------
END MODULE YOE_AERVOLE
