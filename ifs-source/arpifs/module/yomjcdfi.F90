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

MODULE YOMJCDFI

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! RACCSPA*:  Accumulation arrays of spectral variables used
!            in the Digital Filter

REAL(KIND=JPRB), ALLOCATABLE :: RACCSPA3(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: RACCSPA2(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: RACCSPA1(:,:)

! QNORM*: Definition of the norms used for Jc

REAL(KIND=JPRB), ALLOCATABLE :: QNORM1(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: QNORM2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: QNORM3(:,:,:)

! Empirical factors weighting variables in the Jc norm

REAL(KIND=JPRB) :: QNORM_ADJUST_VOR
REAL(KIND=JPRB) :: QNORM_ADJUST_DIV
REAL(KIND=JPRB) :: QNORM_ADJUST_T
REAL(KIND=JPRB) :: QNORM_ADJUST_Q
REAL(KIND=JPRB) :: QNORM_ADJUST_LNSP

INTEGER(KIND=JPIM) :: NSUBDFI
REAL(KIND=JPRB), ALLOCATABLE :: DFICHECK(:)

! RSUMJCDFI*:   Arrays containing the contribution of a given processor
!               to the Jc-DFI cost function (summed in EVCOST)
! RSUMJCDFI1D contribution of the 1D mean wind part of JC

REAL(KIND=JPRB), ALLOCATABLE :: RSUMJCDFI(:)
REAL(KIND=JPRB):: RSUMJCDFI1D

END MODULE YOMJCDFI
