! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE YOMGT1C0

!     Grid point array at time t in the atmospheric physics and dynamics

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE PARDIM1C , ONLY : JPFLEV

IMPLICIT NONE

SAVE

!  === Momentum variables ===

REAL(KIND=JPRB) :: UT0(JPFLEV)
REAL(KIND=JPRB) :: VT0(JPFLEV)

!  === Thermodynamical variables ===

REAL(KIND=JPRB) :: TT0(JPFLEV)
REAL(KIND=JPRB) :: QT0(JPFLEV)

!  === cloud variables ===

REAL(KIND=JPRB) :: WT0(JPFLEV)
REAL(KIND=JPRB) :: ST0(JPFLEV)
REAL(KIND=JPRB) :: AT0(JPFLEV)
REAL(KIND=JPRB) :: RNT0(JPFLEV)
REAL(KIND=JPRB) :: SNT0(JPFLEV)

! === surface pressure ===

REAL(KIND=JPRB) :: SPT0

! === relaxation variables ===

REAL(KIND=JPRB) :: UOBS(JPFLEV)
REAL(KIND=JPRB) :: VOBS(JPFLEV)

REAL(KIND=JPRB) :: TOBS(JPFLEV)
REAL(KIND=JPRB) :: QOBS(JPFLEV)

!     ------------------------------------------------------------------

END MODULE YOMGT1C0
