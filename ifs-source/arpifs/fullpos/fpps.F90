! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPPS(KPROMA,KST,KND,KFLEV,KOPLEV,PIN_GEOPH,POUT_GEOP,PT,&
 & PR,PRESH,PST,PSPPP)

!**** *FPPS*  - FULL-POS pressure post-processing

!     PURPOSE.
!     --------
!        To compute the pressures of a set of surfaces, given their 
!        geopotentials

!**   INTERFACE.
!     ----------
!       *CALL* *FPPS*

!        EXPLICIT ARGUMENTS
!        --------------------

!        * INPUT:
!        KPROMA   : horizontal dimension
!        KST      : start of work
!        KND      : depth of work
!        KFLEV    : number of input levels
!        KOPLEV   : number of output surfaces
!        PIN_GEOPH: input geopotential at half input levels
!        POUT_GEOP: output geopotential
!        PT       : input temperature
!        PR       : input R (moist air constant)
!        PRESH    : input half level hydrostatic pressure
!        PST      : input surface temperature

!        * OUTPUT:
!        PSPPP    : output pressures

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!   ORIGINAL   : 00-08-10 after VEINE
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Aug 2008): move calls to GP.. routines in the caller.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMSTA   , ONLY : RDTDZ1
USE YOMCST   , ONLY : RG

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOPLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN_GEOPH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POUT_GEOP(KPROMA,KOPLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PST(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPPP(KPROMA,KOPLEV) 

!---------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZMDPHIS(KPROMA)         ! Difference (input PHIs - output PHI)
REAL(KIND=JPRB) :: ZPBLO(KPROMA)           ! output Ps for (output orog < input orog)
REAL(KIND=JPRB) :: ZPABV(KPROMA)           ! output Ps for (output orog > input orog)
REAL(KIND=JPRB) :: ZA(KPROMA)              ! =1 if (output orog < input orog) ; else =0
REAL(KIND=JPRB) :: ZTSC(KPROMA)            ! standart temperature for the output surface
REAL(KIND=JPRB) :: ZDTDPHI                 ! standart gradient dT/(g*dZ)

REAL(KIND=JPRB) :: ZSURF_PRES(KPROMA),ZIN_OROG(KPROMA) 

INTEGER(KIND=JPIM) :: IL(KPROMA) ! index of first input half level below output surface
INTEGER(KIND=JPIM) :: JLP, JLEV, JI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------------

#include "pppmer.intfb.h"

!---------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPPS',0,ZHOOK_HANDLE)

!---------------------------------------------------------------------------

!*       1.    COMPUTATION OF SURFACE PRESSURE
!              -------------------------------

ZDTDPHI=RDTDZ1/RG
ZSURF_PRES(KST:KND)=PRESH(KST:KND,KFLEV)
ZIN_OROG(KST:KND)=PIN_GEOPH(KST:KND,KFLEV)

DO JLP=1,KOPLEV

  !*     1.1   Compare surfaces

  ZMDPHIS(KST:KND)=ZIN_OROG(KST:KND)-POUT_GEOP(KST:KND,JLP)
  ZA(KST:KND)=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZMDPHIS(KST:KND)))

  !*     1.2   Case (output orog < input orog) : use msl pressure calculation

  IF (MAXVAL(ZA(KST:KND)) > 0.0_JPRB) THEN
    ZTSC(KST:KND)=PST(KST:KND)-ZDTDPHI*ZMDPHIS(KST:KND)
    CALL PPPMER(KPROMA,KST,KND,ZSURF_PRES,ZMDPHIS,PST,ZTSC,ZPBLO)
  ELSE
    ZPBLO(KST:KND)=0.0_JPRB
  ENDIF

  !*     1.3   Case (output orog > input orog) : P=Po.exp((Phi-Phio)/RT) 

  IF (MINVAL(ZA(KST:KND)) < 1.0_JPRB) THEN

    ! Find nearest model half level below output surface
    IL(KST:KND)=KFLEV
    DO JLEV = KFLEV, 1, -1
      DO JI = KST,KND
        IL(JI) = IL(JI) + &
         & MAX(0.0_JPRB,SIGN(1.0_JPRB,POUT_GEOP(JI,JLP)-PIN_GEOPH(JI,JLEV)))*&
         & (JLEV-IL(JI))  
      ENDDO
    ENDDO

    ! Compute pressure
    DO JI = KST,KND
      ZPABV(JI) = PRESH(JI,IL(JI))*&
       & EXP((PIN_GEOPH(JI,IL(JI))-POUT_GEOP(JI,JLP))/(PR(JI,IL(JI))*PT(JI,IL(JI))))  
    ENDDO

  ELSE

    ZPABV(KST:KND)=0.0_JPRB

  ENDIF

  !*     1.4   Combination 

  PSPPP(KST:KND,JLP)=ZA(KST:KND)*ZPBLO(KST:KND) &
   & +(1.0_JPRB-ZA(KST:KND))*ZPABV(KST:KND)  

ENDDO

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPPS',1,ZHOOK_HANDLE)
END SUBROUTINE FPPS
