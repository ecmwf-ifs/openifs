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

MODULE YEMJK

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! Logical to activate Jk
LOGICAL :: LEJK

! IntegerTruncation of the low resolution LAM geom
INTEGER(KIND=JPIM) :: NSMAXJK

! Coefficients:w to weight the Jk contribution
REAL(KIND=JPRB) :: ALPHAKT ! amplification factor for JK T term
REAL(KIND=JPRB) :: ALPHAKVOR ! amplification factor for JK vorticity term
REAL(KIND=JPRB) :: ALPHAKDIV ! amplification factor for JK divergence term
REAL(KIND=JPRB) :: ALPHAKQ ! amplification factor for JK humidity term
REAL(KIND=JPRB) :: ALPHAKP ! amplification factor for JK Surf. pressure term
REAL(KIND=JPRB) :: PRESINFJK !Pressure level above which JK != 0
REAL(KIND=JPRB) :: PRESUPJK  !Pressure level above which JK has its full amplitude
REAL(KIND=JPRB), ALLOCATABLE :: ZALPHAK(:,:) ! array of amplification factors for JK
                                             !(by level and parameter)


! Norms used for Jk
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: FJKNORM2(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: FJKNORM3(:,:,:)

! Jk Gradient tables
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: VAJKGRAVOR(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: VAJKGRADIV(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: VAJKGRAT(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: VAJKGRAQ(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: VAJKGRASP(:)

! Jk Innovation tables
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: SPJKINVOR(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: SPJKINDIV(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: SPJKINT(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: SPJKINQ(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: SPJKINSP(:)

END MODULE YEMJK
