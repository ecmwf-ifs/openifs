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

MODULE YOMSWE

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Shallow-water equation Williamson et al tests

! ALPHASWE  : Angle ALPHA of the axis of the experiment and the axis of the Earth
! GMUCENSWE : sine  of latitude of the axis of the experiment

REAL(KIND=JPRB) :: ALPHASWE
REAL(KIND=JPRB) :: GMUCENSWE

!     ------------------------------------------------------------------
END MODULE YOMSWE
