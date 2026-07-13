MODULE SOURCE_E_MOD

CONTAINS
SUBROUTINE SOURCE_E(KIDIA,KFDIA,KLON,Z,T,&
 &                  PPHIOC,Z0, WSTAR,XKS,DUDZ,SDS,DTDZ,&
 &                  Q,FST,&
 &                  STOT,DSDE,D_M,&
 &                  D_H,D_E,&
 &                  YDMLM, YDCST)
!
! (C) Copyright 2010- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!***  DETERMINES TKE SOURCE FUNCTION AND DIFFUSION COEFFICIENTS
!
!
!
!
!     TREATMENT OF BUOYANCY EFFECT BY FOLLOWING A
!     FORMULATION  FROM NOH AND KIM (1999) AND SUKORIANSKY ET AL (2005). 
!     THIS IS BASED ON THE TURBULENT RICHARDSON NUMBER RIT=N^2/Q^2, 
!     WITH N THE BRUNT-VAISALA FREQUENCY AND Q THE TURBULENT VELOCITY
!
!
!     AUTHOR: PETER A.E.M. JANSSEN, FEBRUARY 2010
!     ------
!
!
!     PURPOSE.
!     --------
!
!         THE TURBULENT KINETIC ENERGY EQUATION IS SOLVED WITH AN IMPLICIT
!         INTEGRATION IN TIME IN SUB OC_MLM. HERE THE SOURCE FUNCTION 
!         IS DETERMINED, WHERE THE BALANCE INVOLVES TURBULENT DISSIPATION, 
!         SHEAR PRODUCTION, BUOYANCY, LANGMUIR TURBULENT PRODUCTION AND 
!         WAVE BREAKING.   
!
!**   INTERFACE.
!     ----------
!
!         CALL SOURCE_E(KIDIA,KFDIA,KLON,Z,T,PPHIOC,Z0,
!    V                  WSTAR,XKS,DUDZ,SDS,DTDZ,Q,FST,STOT,DSDE,D_M,
!    V                  D_H,D_E,YDMLM)
!         
!         INPUT:
!         -----
!
!         *KIDIA*    INTEGER       FIRST INDEX
!         *KFDIA*    INTEGER       LAST INDEX
!         *KLON*     INTEGER       NUMBER OF GRID POINTS PER PACKET
!         *Z*        REAL          DEPTH
!         *T*        REAL          TEMPERATURE AT DEPTH Z  !!!! in degree C !!!!!
!         *PPHIOC*   REAL          ENERGY FLUX DUE TO OCEAN WAVE DISSIPATION  
!         *Z0*       REAL          ROUGHNESS LENGTH
!         *WSTAR*    REAL          WATER FRICTION VELOCITY
!         *XKS*      REAL          PEAK WAVE NUMBER
!         *DUDZ*     REAL          GRADIENT IN CURRENT
!         *SDS*      REAL          DUDZ . U_STOKES  
!         *DTDZ*     REAL          GRADIENT IN TEMPERATURE
!         *Q*        REAL          VALUE OF TURBULENT VELOCITY
!         
!         OUTPUT:
!         ------
!
!         *FST*      REAL          STOKES-DEPTH PROFILE  
!         *STOT*     REAL          TOTAL SOURCE FUNCTION
!         *DSDE*     REAL          DERIVATIVE OF SOURCE FUNCTION
!         *D_M*      REAL          MOMENTUM DIFFUSION COEFFICIENT
!         *D_H*      REAL          HEAT DIFFUSION COEFFICIENT
!         *D_E*      REAL          TKE DIFFUSION COEFFICIENT
!
!     METHOD.
!     -------
!
!
!     EXTERNALS.
!     ----------
!
!         NO EXTERNALS.
!
!     REFERENCE.
!     ----------
!      
!         OCEAN WAVE EFFECTS ON THE DAILY CYCLE IN SST  BY P.A.E.M. 
!         JANSSEN, 12 OCTOBER 2010.
!
!
!     HEALTH WARNING
!     --------------
!
!         CODE IS WRITTEN ASSUMING DEPTH Z IS POSITIVE (Z => -Z)		
!
!----------------------------------------------------------------------
!
!
USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_MLM,  ONLY : TMLM
USE YOS_CST,  ONLY : TCST

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA, KLON
REAL(KIND=JPRB), INTENT(IN) :: Z
REAL(KIND=JPRB), DIMENSION(KLON), INTENT(IN) :: T, PPHIOC, Z0, WSTAR, XKS
REAL(KIND=JPRB), DIMENSION(KLON), INTENT(IN) :: DUDZ, SDS, DTDZ, Q
REAL(KIND=JPRB), DIMENSION(KLON), INTENT(OUT) :: FST, STOT, DSDE, D_M, D_H, D_E
TYPE(TMLM), INTENT(IN) :: YDMLM
TYPE(TCST), INTENT(IN) :: YDCST

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB), PARAMETER :: ZERO   = 0.0_JPRB
REAL(KIND=JPRB), PARAMETER :: ONE    = 1.0_JPRB
REAL(KIND=JPRB), PARAMETER :: TWO    = 2.0_JPRB
REAL(KIND=JPRB), PARAMETER :: THREE  = 3.0_JPRB

REAL(KIND=JPRB), PARAMETER :: EPS  = 0.00001_JPRB
REAL(KIND=JPRB), PARAMETER :: XNUW = 0.8E-6_JPRB
REAL(KIND=JPRB), PARAMETER :: A_M  = 100._JPRB
REAL(KIND=JPRB), PARAMETER :: A_H  = 80._JPRB 
REAL(KIND=JPRB), PARAMETER :: B_M  = 0.20_JPRB
REAL(KIND=JPRB), PARAMETER :: B_H  = 0.0_JPRB
REAL(KIND=JPRB), PARAMETER :: A_U  = 20._JPRB
REAL(KIND=JPRB), PARAMETER :: PR_1 = 1.2_JPRB
REAL(KIND=JPRB), PARAMETER :: TAU_L = 300._JPRB

REAL(KIND=JPRB) :: ALPHA_W, XN2, XS, XL, RIG, XNSTAR2 
REAL(KIND=JPRB) :: DUSDZ, ZN, FB, DFBDZRHOWM1, RI, XL_M, XL_H
REAL(KIND=JPRB) :: XLM, XLH, BUOY, AA, BB, CC, ZFAC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB), DIMENSION(KLON) :: XL_OBU, SPR, SBR, SLA, SDISS

!--------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SOURCE_E_MOD:SOURCE_E',0,ZHOOK_HANDLE)

ASSOCIATE(RG=>YDCST%RG, B=>YDMLM%B, DML=>YDMLM%DML, F=>YDMLM%F, &
 & R=>YDMLM%R, R1=>YDMLM%R1, R2=>YDMLM%R2, R3=>YDMLM%R3, &
 & RHO=>YDMLM%RHO, SM=>YDMLM%SM, SQ=>YDMLM%SQ, XKAPPA=>YDMLM%XKAPPA, &
 & ZR1=>YDMLM%ZR1, ZR2=>YDMLM%ZR2, ZR3=>YDMLM%ZR3, &
 & ZS=>YDMLM%ZS)

!
!***  1. DETERMINE PARAMETERS.
!     -----------------------
!
DO JL=KIDIA,KFDIA


!!!  ALPHA_W = (T(JL)-273._JPRB)*1.0E-5_JPRB 
  ALPHA_W = T(JL)*1.0E-5_JPRB 
  XN2     = -RG*ALPHA_W*DTDZ(JL)
  XS      = MAX(ABS(DUDZ(JL)),EPS)
!
!***  2. DETERMINE SOURCES.
!     --------------------
! 
  IF (XN2 .GE. ZERO) THEN
     XL = XKAPPA*MIN(TAU_L*WSTAR(JL),(Z+EPS))
  ELSE
     XL = XKAPPA*(Z+EPS)
  ENDIF
           
  FST(JL) = EXP(-TWO*XKS(JL)*Z)
  DUSDZ   = -TWO*XKS(JL)*FST(JL)

  ZN    = Z/Z0(JL)
  FB    = EXP(-ZN)
  DFBDZRHOWM1 = (-ONE/Z0(JL)*FB)/RHO
!
!*    2.1 DETERMINE BUOYANCY PART OF MIXING (NOH & KIM + SUKORIANSKY).
!     ----------------------------------------------------------------
! 
  RI = XN2*XL**2/Q(JL)**2
  IF (XN2 .GE. 0.0_JPRB) THEN
     XL_M = B_M+(ONE-B_M)/SQRT(ONE+A_M*RI)
     XL_H = PR_1*(B_H+(ONE-B_H)/SQRT(ONE+A_H*RI))
  ELSE
     XL_M = ONE-A_U*RI/(ONE-A_U*RI)
     XL_H = PR_1*XL_M
  ENDIF
!
!*    2.2 DETERMINE MIXING LENGTH FOR HEAT AND MOMENTUM.
!     -------------------------------------------------
! 
  XLM = XL*XL_M 
  XLH = XL*XL_H 
!
!
!***  2.3. DETERMINE DIFFUSION COEFFICIENTS.
!     -----------------------------------
!
  D_M(JL)  = XLM*Q(JL)*SM+XNUW
  D_H(JL)  = XLH*Q(JL)*SM+XNUW
  D_E(JL)  = XLM*Q(JL)*SQ+XNUW

  IF (XN2.NE.ZERO) THEN
     XL_OBU(JL) = WSTAR(JL)**3/(XKAPPA*D_H(JL)*XN2)
  ELSE
     XL_OBU(JL) = 1000._JPRB
  ENDIF
!
!*    2.4 DETERMINE LANGMUIR PRODUCTION AND D<PW>/DZ TERM DUE TO BREAKING.
!     -------------------------------------------------------------------
! 
  ZFAC    = SDS(JL)*DUSDZ
  SLA(JL) = D_M(JL)*ZFAC
  CC      = XLM*ZFAC
  SBR(JL) = PPHIOC(JL)*DFBDZRHOWM1
             
!
!*    2.5 DETERMINE SHEAR PRODUCTION AND BUOYANCY TERM.
!     ------------------------------------------------
!
  BUOY = -XLH*XN2
  AA   = SM*(XLM*XS**2+BUOY)

  SPR(JL)   = AA*Q(JL)
  SDISS(JL) = -Q(JL)**3/(B*XLM)
!
!*    2.6 DETERMINE TOTAL SOURCE FUNCTION AND DERIVATIVE.
!     --------------------------------------------------
!
  STOT(JL) = SPR(JL)+SBR(JL)+SLA(JL)+SDISS(JL)
  DSDE(JL) = (AA+CC)/Q(JL)-THREE*Q(JL)/(B*XLM)

ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SOURCE_E_MOD:SOURCE_E',1,ZHOOK_HANDLE)

END SUBROUTINE SOURCE_E

END MODULE SOURCE_E_MOD
