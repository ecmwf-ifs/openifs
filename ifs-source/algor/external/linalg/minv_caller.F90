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

SUBROUTINE MINV_CALLER(LDSCALE,KDIM,PIN,POU)

!**** *MINV_CALLER*   - Inversion of square matrices: interface for MINV

!     Purpose.
!     --------
!      Inversion of a square matrix with dimension KDIM*KDIM: interface for MINV
!      Performs an additional rescaling which is not done in MINV.

!**   Interface.
!     ----------
!        *CALL* *MINV_CALLER(...)

!     Explicit arguments :
!     --------------------
!        LDSCALE          : Activate rescaling                      (IN)
!        KDIM             : Dimension                               (IN)
!        PIN              : Input array                             (IN)
!        POU              : Output array                            (OUT)

!     Implicit arguments :
!     --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation (IDSI)

!     Author.
!     -------
!        K. YESSAD (after SUSI and MINV).
!        Original : July 2017

! Modifications
! -------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL           ,INTENT(IN)     :: LDSCALE
INTEGER(KIND=JPIM),INTENT(IN)     :: KDIM
REAL(KIND=JPRB)   ,INTENT(IN)     :: PIN(KDIM,KDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: POU(KDIM,KDIM)

!     ------------------------------------------------------------------

REAL(KIND=JPRD) :: ZDET, ZEPS
REAL(KIND=JPRB) :: ZDIAG_AVG
REAL(KIND=JPRD) :: ZWORK(2*KDIM)
REAL(KIND=JPRD) :: ZMO(KDIM,KDIM)

INTEGER(KIND=JPIM) :: JL1, JL2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "minv_8.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MINV_CALLER',0,ZHOOK_HANDLE)


IF( LDSCALE )THEN
  ZDIAG_AVG=0.0_JPRB
  DO JL1=1,KDIM
    ZDIAG_AVG = ZDIAG_AVG + ABS(PIN(JL1,JL1))
  ENDDO
  ZDIAG_AVG = ZDIAG_AVG/REAL(KDIM,JPRB)
  ! rescaling
  IF( ZDIAG_AVG > 0.0_JPRB )THEN
    DO JL1=1,KDIM
      DO JL2=1,KDIM
        ZMO(JL1,JL2)=PIN(JL1,JL2)/ZDIAG_AVG
      ENDDO
    ENDDO
  ELSE
    ZMO(:,:)=PIN(:,:)
  ENDIF
ELSE
  ZMO(:,:)=PIN(:,:)
ENDIF

ZEPS=100.0_JPRD*TINY(1.0_JPRD)
CALL MINV_8(ZMO,KDIM,KDIM,ZWORK,ZDET,ZEPS,0,1)

IF( LDSCALE )THEN
  ! rescaling back to original values
  IF( ZDIAG_AVG > 0.0_JPRB )THEN
    DO JL1=1,KDIM
      DO JL2=1,KDIM
        POU(JL1,JL2)=ZMO(JL1,JL2)/ZDIAG_AVG
      ENDDO
    ENDDO
  ELSE
    POU(:,:)=ZMO(:,:)
  ENDIF
ELSE
  POU(:,:)=ZMO(:,:)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MINV_CALLER',1,ZHOOK_HANDLE)
END SUBROUTINE MINV_CALLER
