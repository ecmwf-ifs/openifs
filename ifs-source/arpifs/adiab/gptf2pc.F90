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

SUBROUTINE GPTF2PC(&
 ! --- INPUT ---------------------------------------------------------
 & YDGEOMETRY,YDDYN,YDDYNA,KST,KEN,&
 & PUT9  ,PVT9   ,PTT9  ,PSPDT9  ,PSVDT9  ,&
 & PSPT9 ,PNHXT9 ,&
 & PUT0  ,PVT0   ,PTT0  ,PSPDT0  ,PSVDT0  ,&
 & PSPT0 ,PNHXT0    ,&
 ! --- INPUT/OUTPUT --------------------------------------------------
 & PCUT9 ,PCVT9  ,PCTT9 ,PCSPDT9 ,PCSVDT9 ,&
 & PCSPT9,PCNHXT9)

!**** *GPTF2PC* - Timefilter part 2

!     Purpose.
!     --------
!           Performs part 2 of the time-filtering for the t-dt array.
!           Routine called in 3TL PC schemes.
            
!**   Interface.
!     ----------
!        *CALL* *GPTF2PC(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KST                    : start of work
!         KEN                    : depth of work
!         P(U,V,T,SPD,SVD,SP)T9  : pre-filtered quantities at time t-dt
!         PNHXT9                 : pre-filtered Xterm at time t-dt
!         P(U,V,T,SPD,SVD,SP)T0  : quantities at time t
!         PNHXT0                 : Xterm at time t

!        INPUT/OUTPUT:
!         PC(U,V,T,SPD,SVD,SP)T9 : quantities at time t-dt
!         PCNHXT9                : Xterm at time t-dt

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        LACE internal documentation

!     Author.
!     -------
!        Jozef Vivoda (based on GPTF1)  *ECMWF*
!        Original : 21-Feb-2005

! Modifications
! -------------
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (June 2017): Introduce NHQE model.
! End Modifications
!------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT3       , ONLY : NSTEP
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSVDT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSVDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCUT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCVT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSPDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSVDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSPT9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL,JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF2PC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, REPS2=>YDDYN%REPS2)

!     ------------------------------------------------------------------

!*       1. PERFORM TIME FILTER (PART 1)
!           ----------------------------

IF (YDDYNA%LTWOTL.OR.(.NOT.YDDYNA%LPC_FULL).OR.(NCURRENT_ITER/=0)) THEN
  CALL ABOR1(' GPTF2PC: Case where GPTF2PC should not be called!')
ENDIF

IF( NSTEP <= 0 )THEN

  ! U,V,T,log(Pi_s)
  DO JL=1,NFLEVG
    DO JK=KST,KEN
      PCUT9(JK,JL)=PUT0(JK,JL)
      PCVT9(JK,JL)=PVT0(JK,JL)
      PCTT9(JK,JL)=PTT0(JK,JL)
    ENDDO
  ENDDO
  DO JK=KST,KEN
    PCSPT9(JK)=PSPT0(JK)
  ENDDO

  ! NH variables
  IF(YDDYNA%LNHEE)THEN
    DO JL=1,NFLEVG
      DO JK=KST,KEN
        PCSPDT9(JK,JL)=PSPDT0(JK,JL)
      ENDDO
    ENDDO
  ENDIF
  IF(YDDYNA%LNHDYN)THEN
    DO JL=1,NFLEVG
      DO JK=KST,KEN
        PCSVDT9(JK,JL)=PSVDT0(JK,JL)
      ENDDO
    ENDDO
  ENDIF
  IF(YDDYNA%LNHDYN.AND.YDDYNA%NVDVAR==4)THEN
    DO JL=1,NFLEVG
      DO JK=KST,KEN
        PCNHXT9(JK,JL)=PNHXT0(JK,JL)
      ENDDO
    ENDDO
  ENDIF

ELSE

  ! U,V,T,log(Pi_s)
  DO JL=1,NFLEVG
    DO JK=KST,KEN
      PCUT9(JK,JL)=PUT9(JK,JL)+REPS2*PUT0(JK,JL)
      PCVT9(JK,JL)=PVT9(JK,JL)+REPS2*PVT0(JK,JL)
      PCTT9(JK,JL)=PTT9(JK,JL)+REPS2*PTT0(JK,JL)
    ENDDO
  ENDDO
  DO JK=KST,KEN
    PCSPT9(JK)=PSPT9(JK)+REPS2*PSPT0(JK)
  ENDDO

  ! NH variables
  IF(YDDYNA%LNHEE)THEN
    DO JL=1,NFLEVG
      DO JK=KST,KEN
        PCSPDT9(JK,JL)=PSPDT9(JK,JL)+REPS2*PSPDT0(JK,JL)
      ENDDO
    ENDDO
  ENDIF
  IF(YDDYNA%LNHDYN)THEN
    DO JL=1,NFLEVG
      DO JK=KST,KEN
        PCSVDT9(JK,JL)=PSVDT9(JK,JL)+REPS2*PSVDT0(JK,JL)
      ENDDO
    ENDDO
  ENDIF
  IF(YDDYNA%LNHDYN.AND.YDDYNA%NVDVAR==4)THEN
    DO JL=1,NFLEVG
      DO JK=KST,KEN
        PCNHXT9(JK,JL)=PNHXT9(JK,JL)+REPS2*PNHXT0(JK,JL)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF2PC',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF2PC
