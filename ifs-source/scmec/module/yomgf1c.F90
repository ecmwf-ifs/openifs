! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE YOMGF1C

!*    grid point array for large scale forcing of the one column model

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE PARDIM1C , ONLY : JPFLEV

IMPLICIT NONE

SAVE

!  === Geostrophic Wind ===

REAL(KIND=JPRB) :: UG0(JPFLEV)
REAL(KIND=JPRB) :: VG0(JPFLEV)

!  === Vertical Velocity ===

REAL(KIND=JPRB), target :: VVEL0(1,JPFLEV)

!  === horizontal advection ===

REAL(KIND=JPRB) :: UADV(JPFLEV)
REAL(KIND=JPRB) :: VADV(JPFLEV)
REAL(KIND=JPRB) :: TADV(JPFLEV)
REAL(KIND=JPRB) :: QADV(JPFLEV)

REAL(KIND=JPRB) :: etadotdpdeta(0:JPFLEV)

!     ------------------------------------------------------------------

END MODULE YOMGF1C
