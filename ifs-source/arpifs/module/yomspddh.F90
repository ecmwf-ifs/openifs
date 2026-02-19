! (C) Copyright 1989- Meteo-France.

MODULE YOMSPDDH

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Spectral arrays used in Horizontal Diffusion (HD):
!     and Semi-Implicit (SI) DDH tendencies.

! -   ARRAYS

TYPE :: TSPDDH
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_VOR(:,:) ! for vorticity
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_DIV(:,:) ! for divergence
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_T(:,:)   ! for temperature
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_SPD(:,:) ! for NH pressure departure
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_SVD(:,:) ! for NH vertical divergence
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_SNHX(:,:)! for NH "X" part of divergence
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDHD_GFL(:,:,:) ! for specific humidity

REAL(KIND=JPRB),ALLOCATABLE:: SPTNDSI_VOR(:,:) ! for vorticity
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDSI_DIV(:,:) ! for divergence
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDSI_T(:,:)   ! for temperature
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDSI_SPD(:,:) ! for NH pressure departure
REAL(KIND=JPRB),ALLOCATABLE:: SPTNDSI_SVD(:,:) ! for NH vertical divergence
END TYPE TSPDDH

!!TYPE(TSPDDH), POINTER :: YRSPDDH => NULL()

!     ------------------------------------------------------------------
END MODULE YOMSPDDH
