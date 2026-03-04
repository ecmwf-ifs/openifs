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

SUBROUTINE SULEGA(YDCSGLEG,YDDIM)

!**** *SULEGA * - initialize the Legendre polynomials

!     Purpose.
!     --------
!           Initialize COMMON YOMLEG

!**   Interface.
!     ----------
!        *CALL* *SULEGA*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!              COMMON YOMLEG

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15 (originally named SULEG)

!     Modifications.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE YOMLEG   , ONLY : TCSGLEG
USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI, RA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCSGLEG), INTENT(INOUT) :: YDCSGLEG
TYPE(TDIM)   , INTENT(IN) :: YDDIM

INTEGER(KIND=JPIM) :: JGL
REAL(KIND=JPRD) :: ZZD, ZZDS
REAL(KIND=JPRD), ALLOCATABLE :: ZMU(:), ZPHI(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "trans_inq.h"  

!     ------------------------------------------------------------------ 
IF (LHOOK) CALL DR_HOOK('SULEGA',0,ZHOOK_HANDLE)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG, &
 & NDGSUR=>YDDIM%NDGSUR, NRESOL=>YDDIM%NRESOL)
!     ------------------------------------------------------------------ 

!*       1.    Allocate space (moved from SUALLO)
!              ----------------------------------

! ky: allocated if LELAM=T too to avoid problems but they should not be used under LELAM=T.
!     RMU and RSQM2 appear in SLCSET dummy arg interface, but are used under SLCSET for LELAM=F only.
!     RLATIG is used in SUGRIB for LELAM=T too.

ALLOCATE(YDCSGLEG%RW    (NDGLG))
ALLOCATE(ZMU(NDGLG))
ALLOCATE(ZPHI(NDGLG))
ALLOCATE(YDCSGLEG%RMU   (NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%R1MU2 (NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%R1MUI (NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%R1MUA (NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%RSQM2 (NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%R1QM2 (NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%RACTHE(NDGSAG:NDGENG))
ALLOCATE(YDCSGLEG%RLATIG(NDGSAG:NDGENG))

!*       2.    Gaussian latitudes and weights 
!              ------------------------------ 

CALL TRANS_INQ(KRESOL=NRESOL,PMU=ZMU(1:NDGLG),PGW=YDCSGLEG%RW)
YDCSGLEG%RMU(1:NDGLG) = ZMU(1:NDGLG)

!*       3.    Computes related arrays
!              -----------------------

DO JGL=1,NDGLG
  ZZD = 1.0_JPRD - ZMU(JGL)*ZMU(JGL)
  ZZDS = SQRT(ZZD)
  ZPHI(JGL) = ASIN(ZMU(JGL))
  YDCSGLEG%R1MU2(JGL) = ZZD
  YDCSGLEG%R1MUI(JGL) = REAL(1.0_JPRD/ZZD,JPRB)
  YDCSGLEG%R1MUA(JGL) = REAL(1.0_JPRD/(REAL(RA,JPRD)*ZZD),JPRB)
  YDCSGLEG%RSQM2(JGL) = ZZDS
  YDCSGLEG%R1QM2(JGL) = REAL(1.0_JPRD/ZZDS,JPRB)
  YDCSGLEG%RACTHE(JGL) = REAL(1.0_JPRD/(REAL(RA,JPRD)*ZZDS),JPRB)
  YDCSGLEG%RLATIG(JGL) = ZPHI(JGL)
ENDDO
DEALLOCATE(ZMU)

!*       4.    Computes poles values.
!              ---------------------

IF(NDGSUR >= 1)THEN

!DIR$ IVDEP
!OCL NOVREC
  DO JGL=1,NDGSUR
    YDCSGLEG%RMU(1-JGL)=YDCSGLEG%RMU(JGL)
    YDCSGLEG%R1MU2(1-JGL)=YDCSGLEG%R1MU2(JGL)
    YDCSGLEG%R1MUI(1-JGL)=YDCSGLEG%R1MUI(JGL)
    YDCSGLEG%R1MUA(1-JGL)=YDCSGLEG%R1MUA(JGL)
    YDCSGLEG%RSQM2(1-JGL)=YDCSGLEG%RSQM2(JGL)
    YDCSGLEG%R1QM2(1-JGL)=YDCSGLEG%R1QM2(JGL)
    YDCSGLEG%RACTHE(1-JGL)=YDCSGLEG%RACTHE(JGL)
    YDCSGLEG%RLATIG(1-JGL)=REAL(RPI,JPRD)-ZPHI(JGL)
  ENDDO

!DIR$ IVDEP
!OCL NOVREC
  DO JGL=1,NDGSUR
    YDCSGLEG%RMU(NDGLG+JGL)=YDCSGLEG%RMU(NDGLG+1-JGL)
    YDCSGLEG%R1MU2(NDGLG+JGL)=YDCSGLEG%R1MU2(NDGLG+1-JGL)
    YDCSGLEG%R1MUI(NDGLG+JGL)=YDCSGLEG%R1MUI(NDGLG+1-JGL)
    YDCSGLEG%R1MUA(NDGLG+JGL)=YDCSGLEG%R1MUA(NDGLG+1-JGL)
    YDCSGLEG%RSQM2(NDGLG+JGL)=YDCSGLEG%RSQM2(NDGLG+1-JGL)
    YDCSGLEG%R1QM2(NDGLG+JGL)=YDCSGLEG%R1QM2(NDGLG+1-JGL)
    YDCSGLEG%RACTHE(NDGLG+JGL)=YDCSGLEG%RACTHE(NDGLG+1-JGL)
    YDCSGLEG%RLATIG(NDGLG+JGL)=-REAL(RPI,JPRD)-ZPHI(NDGLG+1-JGL)
  ENDDO

ENDIF
DEALLOCATE(ZPHI)
!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SULEGA',1,ZHOOK_HANDLE)
END SUBROUTINE SULEGA
