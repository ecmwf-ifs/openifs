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

MODULE TYPE_FPOSBUF

USE PARKIND1  ,ONLY : JPRB
USE YOMFP4L , ONLY : TRQFP

IMPLICIT NONE

SAVE

! Post-prpcessing data derived type :
! =================================

! YRQPHY: fields request
! FPBUF : cache-blocked data array

TYPE FPOSBUF


TYPE(TRQFP)               :: YRQPHY
REAL(KIND=JPRB)    , ALLOCATABLE :: FPBUF(:,:,:)

END TYPE FPOSBUF

END MODULE TYPE_FPOSBUF
