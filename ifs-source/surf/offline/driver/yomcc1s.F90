MODULE YOMCC1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!*    Values for seasonally varying fields

REAL(KIND=JPRB),ALLOCATABLE,TARGET:: GPCC(:,:,:)

! -   Splitting of the array GPCC

!     VCALB*  -  ALBEDO  
!     VCLAI*  -  LAI 
!     VCVEG*  -  VEGETATION COVERAGE
 

REAL(KIND=JPRB),POINTER:: VCALB(:,:)
REAL(KIND=JPRB),POINTER:: VCLAIL(:,:)
REAL(KIND=JPRB),POINTER:: VCLAIH(:,:)
REAL(KIND=JPRB),POINTER:: VCFWET(:,:)
REAL(KIND=JPRB),POINTER:: VCVEG(:,:)
REAL(KIND=JPRB),POINTER:: VCLAI(:,:)

END MODULE YOMCC1S
