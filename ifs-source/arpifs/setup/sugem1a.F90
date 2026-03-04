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

SUBROUTINE SUGEM1A(YDGEOMETRY,KSUPERSEDE)

!**** *SUGEM1A*  - Initialize geometry parameters - first part

!     Purpose.
!     --------
!           Initialize geometry: compute NMENG and NDGLU.

!**   Interface.
!     ----------
!        *CALL* *SUGEM1A

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        See #include below

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
!      Original : 87-10-15

!     Modifications.
!     --------------
!      K. Yessad (Sept 2008): Gaussian lat and weights => dummy arg
!      K. Yessad (Sept 2008): Prune conf 951.
!      K. Yessad (Nov 2008): rename arp/SUGAW into arp/SUGAWA.
!      K. Yessad (Aug 2009): prune conf 912, externalise conf 911.
!      K. Yessad (May 2012): simplifications; namelist variables in SUGEM_NAML.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      K. Yessad (Dec 2016): Prune obsolete options.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMARG   , ONLY : NUMEN, NSUPERSEDE
USE YOMCT0   , ONLY : LALLOPR, LECMWF
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSUPERSEDE
INTEGER(KIND=JPIM) :: IDGLU(0:YDGEOMETRY%YRDIM%NSMAX,YDGEOMETRY%YRDIM%NDGSAG:YDGEOMETRY%YRDIM%NDGNH)
INTEGER(KIND=JPIM) :: IDIF(1:(YDGEOMETRY%YRDIM%NDGLG+1)/2)
INTEGER(KIND=JPIM) :: IISMAX,JGL,JM
LOGICAL :: LL_FROM_CL, LL_RECOMPUTE_NMENG, LLP, LLSUPERSEDE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

#include "trans_inq.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGEM1A',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGNH=>YDDIM%NDGNH, &
 & NDGSAG=>YDDIM%NDGSAG, NDGSUR=>YDDIM%NDGSUR, NRESOL=>YDDIM%NRESOL, &
 & NSMAX=>YDDIM%NSMAX, &
 & RNLGINC=>YDGEM%RNLGINC)
!     ------------------------------------------------------------------

!*       0.    Allocations (moved from SUALLO)
!              -------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR

!! LELAM code for SUEGEM1A: IISMAX = NMSMAX
IISMAX = NSMAX

ALLOCATE(YDGEM%NMEN   (NDGSAG:NDGENG))
IF(LLP)WRITE(NULOUT,9) 'NMEN     ',SIZE(YDGEM%NMEN     ),SHAPE(YDGEM%NMEN     )
ALLOCATE(YDGEM%NMENG  (NDGSAG:NDGENG))
IF(LLP)WRITE(NULOUT,9) 'NMENG    ',SIZE(YDGEM%NMENG    ),SHAPE(YDGEM%NMENG    )
ALLOCATE(YDGEM%NDGLU  (0:IISMAX))
IF(LLP)WRITE(NULOUT,9) 'NDGLU    ',SIZE(YDGEM%NDGLU    ),SHAPE(YDGEM%NDGLU    )

!     ------------------------------------------------------------------

!*       1.    Get NMENG and NDGLU from transform setup.
!              -----------------------------------------

CALL TRANS_INQ(KRESOL=NRESOL,KNMENG=YDGEM%NMENG(1:NDGLG),KDGLU=YDGEM%NDGLU(0:NSMAX))
IF (NDGSUR > 0) THEN
!DIR$ IVDEP
  DO JGL=1,NDGSUR
    YDGEM%NMENG(1-JGL)=YDGEM%NMENG(JGL)
    YDGEM%NMENG(NDGLG+JGL)=YDGEM%NMENG(NDGLG+1-JGL)
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    Correct NMENG and NDGLU from file frame if required.
!              ------------------------------------------------------

IF (PRESENT(KSUPERSEDE)) THEN
  LLSUPERSEDE=(KSUPERSEDE == 1)
ELSE
  LLSUPERSEDE=(NSUPERSEDE == 1)
ENDIF

LL_FROM_CL=(.NOT.LECMWF .AND. LLSUPERSEDE)

IF (LL_FROM_CL) THEN

  ! NMENG taken from file frame:

  ! First we control that the value of RNLGINC correspond to what is in file
  ! if not we overwrite NMENG with the input file values
  ! but we signal this mismatch by setting RNLGINC=999.

  DO JGL = 1, (NDGLG+1)/2
    IDIF(JGL)=YDGEM%NMENG(JGL)-NUMEN(JGL)
  ENDDO
  LL_RECOMPUTE_NMENG=.FALSE.
  IF (ANY(IDIF(1:(NDGLG+1)/2) /= 0)) THEN
    WRITE(NULOUT,'('' SUGEM1A: NMENG IN FILE AND RNLGINC IN NAMELIST OR FILE FRAME MISMATCH'')')  
    RNLGINC=999._JPRB
    WRITE(NULOUT,FMT='('' RNLGINC RESET TO '',E13.7, '' IN SUGEM1A '')') RNLGINC
    LL_RECOMPUTE_NMENG=.TRUE.
  ENDIF

  IF (LL_RECOMPUTE_NMENG) THEN
    ! NMENG got from TRANS_INQ has not the right value; copy NUMEN on NMENG.
    YDGEM%NMENG(1:NDGLG)=NUMEN(1:NDGLG)
    IF (NDGSUR > 0) THEN
!DIR$ IVDEP
      DO JGL=1,NDGSUR
        YDGEM%NMENG(1-JGL)=YDGEM%NMENG(JGL)
        YDGEM%NMENG(NDGLG+JGL)=YDGEM%NMENG(NDGLG+1-JGL)
      ENDDO
    ENDIF
  ENDIF

  ! Compute NDGLU from NMENG:
  DO JGL=NDGSAG,NDGNH
    DO JM=0,YDGEM%NMENG(JGL)
      IDGLU(JM,JGL)=1
    ENDDO
    DO JM=YDGEM%NMENG(JGL)+1,NSMAX
      IDGLU(JM,JGL)=0
    ENDDO
  ENDDO
  DO JM=0,NSMAX
    YDGEM%NDGLU(JM)=0
    DO JGL=1,NDGNH
      YDGEM%NDGLU(JM)=YDGEM%NDGLU(JM)+IDGLU(JM,JGL)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*       3.    Printings.
!              ----------

IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(NULOUT,*) ' '
  WRITE(NULOUT,*) ' --- Printings in SUGEM1A:'
  WRITE(NULOUT,FMT='('' (JGL,NMENG) '')')
  WRITE(NULOUT,FMT='(8(1X,''('',I4,I5,'')''))') (JGL,YDGEM%NMENG(JGL),JGL=NDGSAG,NDGENG)  
  WRITE(NULOUT,FMT='('' (JM,NDGLU) '')')
  WRITE(NULOUT,FMT='(10(1X,''('',I4,I4,'')''))') (JM,YDGEM%NDGLU(JM),JM=0,NSMAX)  
ENDIF

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGEM1A',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEM1A
