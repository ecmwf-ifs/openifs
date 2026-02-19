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

MODULE YOMDYNA_STATIC

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE
! ----------------------------------------------------------------------



! ------ Other diffusive processes ------------------------------------------

! LGRADGP    : use grid-point hor. derivatives

LOGICAL :: LGRADGP

! LDRY_ECMWF  : .TRUE. = compute grad(RT) in TL/AD without moist increments

LOGICAL :: LDRY_ECMWF

CONTAINS
SUBROUTINE SUDYNA_STATIC
USE YOMLUN       , ONLY : NULNAM,NULOUT
!------------------------------------------------------------------

IMPLICIT NONE


NAMELIST /NAMDYNA_STATIC/ LGRADGP,LDRY_ECMWF
#include "posnam.intfb.h"

LGRADGP=.FALSE.
LDRY_ECMWF=.FALSE.

CALL POSNAM(NULNAM,'NAMDYNA_STATIC')
READ(NULNAM,NAMDYNA_STATIC)



WRITE(UNIT=NULOUT,FMT='('' Printings of YOMDYNA_STATIC variables '')')
WRITE(UNIT=NULOUT,FMT='('' LGRADGP = '',L2)') LGRADGP
WRITE(UNIT=NULOUT,FMT='('' LDRY_ECMWF= '',L2)') LDRY_ECMWF
END SUBROUTINE SUDYNA_STATIC

! ----------------------------------------------------------------------
END MODULE YOMDYNA_STATIC
