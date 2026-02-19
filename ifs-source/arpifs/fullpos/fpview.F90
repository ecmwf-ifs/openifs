! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPVIEW(KPROMA,KST,KND,KOPLEV,PSP1,PSP2,PRES2,PDPBLC,PALFA)

!**** *FPVIEW* - FULLPOS Vertical Interpolation Entry Weights.

!     PURPOSE.
!     --------
!           Compute weights for vertical interpolations combining 2 profiles

!**   INTERFACE.
!     ----------
!        *CALL* *FPVIEW*

!        EXPLICIT ARGUMENTS :
!        --------------------
!        KPROMA : horizontal dimension
!        KST    : starting value of work.
!        KND    : depth of work.
!        KOPLEV : number of output levels.
!        PSP1   : surface pressure 1
!        PSP2   : surface pressure 2 
!        PDPBLC : critical thickness of PBL (in Pa)
!        PRES2  : pressures 2.
!        PALFA  : Weights for vertical profiles.

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE* 

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 00-08-10 After VEINE
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOPLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP2(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES2(KPROMA,KOPLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDPBLC 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALFA(KPROMA,KOPLEV) 
REAL(KIND=JPRB) :: ZALT(KPROMA)  ! =1 if (|Ps1 -Ps2| > PDPBLC) ; else =0
REAL(KIND=JPRB) :: ZPBOT(KPROMA) ! Pressure on bottom of PBL
REAL(KIND=JPRB) :: ZPTOP(KPROMA) ! Pressure on top of PBL

INTEGER(KIND=JPIM) :: JLEV, JI

REAL(KIND=JPRB) :: ZARG, ZXX, ZEPS

REAL(KIND=JPRB) :: ZFXX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ------------------------------------------------------------------
! * DEFINE REAL FUNCTION FOR CALCULATION OF WEIGHTS ALFA

ZFXX(ZARG) = 3._JPRB*ZARG*ZARG - 2._JPRB*ZARG*ZARG*ZARG

!*       3.    WEIGHTS ALFA FOR 2 VERTICAL PROFILES.
!              -------------------------------------

IF (LHOOK) CALL DR_HOOK('FPVIEW',0,ZHOOK_HANDLE)
DO JI=KST,KND
  ZALT(JI) = MAX(0.0_JPRB,SIGN(1.0_JPRB,ABS(PSP1(JI)-PSP2(JI))-PDPBLC))
  ZPBOT(JI) = ZALT(JI)*MIN(PSP1(JI),PSP2(JI)) &
   & +(1.0_JPRB-ZALT(JI))*(MAX(PSP1(JI),PSP2(JI))-PDPBLC)  
  ZPTOP(JI) = ZALT(JI)*(ZPBOT(JI)-PDPBLC)+(1-ZALT(JI)) *&
   & (MIN(PSP1(JI),PSP2(JI))-PDPBLC)  
ENDDO

ZEPS=0.001_JPRB
DO JLEV=1,KOPLEV
  DO JI=KST,KND
    ZXX = (ZPBOT(JI)-PRES2(JI,JLEV))/MAX((ZPBOT(JI)-ZPTOP(JI)),ZEPS)
    PALFA(JI,JLEV) = MAX(0.0_JPRB,SIGN(1.0_JPRB,ZPBOT(JI)-PRES2(JI,JLEV))) *&
     & (1.0_JPRB + MAX(0.0_JPRB,SIGN(1.0_JPRB,PRES2(JI,JLEV)-ZPTOP(JI)))*&
     & (ZFXX(ZXX) - 1.0_JPRB))  
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('FPVIEW',1,ZHOOK_HANDLE)
END SUBROUTINE FPVIEW
