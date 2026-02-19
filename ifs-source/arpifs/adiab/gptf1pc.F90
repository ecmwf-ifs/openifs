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

SUBROUTINE GPTF1PC(&
 ! --- INPUT ---------------------------------------------------------
 & YDGEOMETRY,YDDYN,YDDYNA,KST,KEN,&
 & PCUT9 ,PCVT9  ,PCTT9 ,PCSPDT9 ,PCSVDT9 ,&
 & PCSPT9,PCNHXT9,&
 & PUT0  ,PVT0   ,PTT0  ,PSPDT0  ,PSVDT0  ,&
 & PSPT0 ,PNHXT0    ,&
 ! --- INPUT/OUTPUT --------------------------------------------------
 & PUT9  ,PVT9   ,PTT9  ,PSPDT9  ,PSVDT9  ,&
 & PSPT9 ,PNHXT9)

!**** *GPTF1PC* - Timefilter part 1

!     Purpose.
!     --------
!           Performs part 1 of the time-filtering for the t-dt array.
!           Routine called in 3TL PC schemes.
           
!           - leap-frog: px9=eps1*pcx9+(1-eps1-eps2)*px0
!           - sl2tl    : no action

!**   Interface.
!     ----------
!        *CALL* *GPTF1PC(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KST                    : start of work
!         KEN                    : depth of work
!         PC(U,V,T,SPD,SVD,SP)T9 : quantities at time t-dt
!         P(U,V,T,SPD,SVD,SP)T0  : quantities at time t                        
!         PCNHXT9                : Xterm at time t-dt
!         PNHXT0                 : Xterm at time t

!        INPUT/OUTPUT:
!         P(U,V,T,SPD,SVD,SP)T9  : pre-filtered quantities at time t
!         PNHXT9                 : pre-filtered Xterm at time t

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
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCTT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSPDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSVDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSPT9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSVDT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSVDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPT9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZREST
INTEGER(KIND=JPIM) :: JL,JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF1PC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, REPS1=>YDDYN%REPS1, REPS2=>YDDYN%REPS2)
!     ------------------------------------------------------------------

!*       1. PERFORM TIME FILTER (PART 1)
!           ----------------------------

IF (YDDYNA%LTWOTL.OR.(.NOT.YDDYNA%LPC_FULL).OR.(NCURRENT_ITER/=0)) THEN
  CALL ABOR1(' GPTF1PC: Case where GPTF1PC should not be called!')
ENDIF

ZREST = 1.0_JPRB-(REPS1+REPS2)

! U,V,T,log(Pi_s)
DO JL=1,NFLEVG
  DO JK=KST,KEN
    PUT9(JK,JL)=REPS1*PCUT9(JK,JL)+ZREST*PUT0(JK,JL)
    PVT9(JK,JL)=REPS1*PCVT9(JK,JL)+ZREST*PVT0(JK,JL)
    PTT9(JK,JL)=REPS1*PCTT9(JK,JL)+ZREST*PTT0(JK,JL)
  ENDDO
ENDDO
DO JK=KST,KEN
  PSPT9(JK) = REPS1*PCSPT9(JK)+ZREST*PSPT0(JK)
ENDDO

! NH variables
IF(YDDYNA%LNHEE)THEN
  DO JL=1,NFLEVG
    DO JK=KST,KEN
      PSPDT9(JK,JL)=REPS1*PCSPDT9(JK,JL)+ZREST*PSPDT0(JK,JL)
    ENDDO
  ENDDO
ENDIF
IF(YDDYNA%LNHDYN)THEN
  DO JL=1,NFLEVG
    DO JK=KST,KEN
      PSVDT9(JK,JL)=REPS1*PCSVDT9(JK,JL)+ZREST*PSVDT0(JK,JL)
    ENDDO
  ENDDO
ENDIF
IF(YDDYNA%LNHDYN.AND.YDDYNA%NVDVAR==4)THEN
  DO JL=1,NFLEVG
    DO JK=KST,KEN
      PNHXT9(JK,JL)=REPS1*PCNHXT9(JK,JL)+ZREST*PNHXT0(JK,JL)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF1PC',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF1PC
