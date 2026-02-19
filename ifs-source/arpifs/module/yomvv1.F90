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

MODULE YOMVV1

USE PARKIND1, ONLY : JPRB
USE PARKIND2, ONLY : JPRH
USE PARDIM,   ONLY : JPMXLE

IMPLICIT NONE

SAVE
!-----------------------------------------------------------------------------

! * Variables used to define the vertical hybrid coordinate:
!   On the half layer "lbar":
!   prehyd(lbar)=DVALH(lbar)+DVBH(lbar)*prehyds
!   (prehyd is the hydrostatic pressure on the half layer "lbar",
!   prehyds is the surface hydrostatic pressure).
!   VVP00 is a reference surface pressure.
!   (DVALH/VVP00) and DVBH vary between 0 and 1.

REAL(KIND=JPRB) :: VVP00
REAL(KIND=JPRH) :: DVALH(0:JPMXLE)
REAL(KIND=JPRH) :: DVBH(0:JPMXLE)

!-----------------------------------------------------------------------------
END MODULE YOMVV1
