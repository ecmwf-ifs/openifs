! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE YOMGT1C1

!     Grid point array at time t+dt in the atmospheric physics and dynamics

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE PARDIM1C , ONLY : JPFLEV

IMPLICIT NONE

SAVE

!  === Momentum variables ===

REAL(KIND=JPRB) :: UT1(JPFLEV)
REAL(KIND=JPRB) :: VT1(JPFLEV)

!  === Thermodynamical variables ===

REAL(KIND=JPRB) :: TT1(JPFLEV)
REAL(KIND=JPRB) :: QT1(JPFLEV)

!  === cloud variables ===

REAL(KIND=JPRB) :: WT1(JPFLEV)
REAL(KIND=JPRB) :: ST1(JPFLEV)
REAL(KIND=JPRB) :: AT1(JPFLEV)
REAL(KIND=JPRB) :: RNT1(JPFLEV)
REAL(KIND=JPRB) :: SNT1(JPFLEV)

! === surface pressure ===

REAL(KIND=JPRB) :: SPT1

!     ------------------------------------------------------------------

END MODULE YOMGT1C1
