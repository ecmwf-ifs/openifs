! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FFT
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE BLUESTEIN_MOD ,ONLY : FFTB_TYPE

! Module for Fourier transforms.

IMPLICIT NONE

SAVE

TYPE FFT_TYPE
  REAL(KIND=JPRB)   ,ALLOCATABLE :: TRIGS(:,:) ! list of trigonometric function values
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFAX(:,:)  ! list of factors of truncation
  LOGICAL                        :: LBLUESTEIN=.FALSE. ! logical indicating whether some
                                                       ! latitudes require bluestein algorithm
                                                       ! with prime factors that are not 2,3,or 5
  LOGICAL,ALLOCATABLE            :: LUSEFFT992(:) ! describes which FFT algorithm to be used
                                                  ! T=use FFT992 F=use bluestein
END TYPE FFT_TYPE

TYPE(FFT_TYPE),ALLOCATABLE,TARGET :: FFT_RESOL(:)
TYPE(FFT_TYPE),POINTER     :: T

TYPE(FFTB_TYPE),ALLOCATABLE,TARGET :: FFTB_RESOL(:)
TYPE(FFTB_TYPE),POINTER    :: TB

END MODULE TPM_FFT
