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

MODULE YOELWCONST

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!*
!      /YOELWCONST/ - CONSTANT USED BY LW RADIATION

REAL(KIND=JPRB) :: RCH4A
REAL(KIND=JPRB) :: RCH4B
REAL(KIND=JPRB) :: RCN2OA
REAL(KIND=JPRB) :: RCN2OB

!     ------------------------------------------------------------------

END MODULE YOELWCONST
