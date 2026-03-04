! (C) Copyright 1989- Meteo-France.

MODULE YOMGPDDH

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    GP Buffer used in Horizontal Diffusion (HD):
!     and Semi-Implicit (SI) DDH tendencies.

! -   ARRAYS

TYPE :: TGPDDH
REAL(KIND=JPRB),ALLOCATABLE:: GMVTNDHD_DDH(:,:,:,:) !
REAL(KIND=JPRB),ALLOCATABLE:: GFLTNDHD_DDH(:,:,:,:) !
REAL(KIND=JPRB),ALLOCATABLE:: GMVTNDSI_DDH(:,:,:,:) !
REAL(KIND=JPRB),ALLOCATABLE:: GMVTNDSL_DDH(:,:,:,:) !
REAL(KIND=JPRB),ALLOCATABLE:: GFLTNDSL_DDH(:,:,:,:) !
END TYPE TGPDDH

!!TYPE(TGPDDH), POINTER :: YRGPDDH => NULL()

!     ------------------------------------------------------------------
END MODULE YOMGPDDH
