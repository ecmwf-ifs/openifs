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

MODULE YOMVSPLIP

USE PARKIND1 , ONLY : JPRB
IMPLICIT NONE

SAVE

! * Variables related to use of spline cubic vertical interpolations.
! RVSPTRI,RVSPC      : are used to re-profile the field to be interpolated
!                      in routine VSPLTRANS.
! RFVV: is used in the computation of the vertical weights.

TYPE TVSPLIP
REAL(KIND=JPRB),ALLOCATABLE :: RVSPTRI(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RVSPC(:)
REAL(KIND=JPRB),ALLOCATABLE :: RFVV(:,:,:)
END TYPE TVSPLIP

END MODULE YOMVSPLIP
