! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLOUDSC2 (YDEPHLI,YDPHNC,YDECLD,YDECLDP , &
 & KIDIA , KFDIA , KLON , KTDIA , KLEV , LDRAIN1D, &
 & PTSPHY, &
 & PAPHP1 , PAPP1 ,  PQM1  , PQS  , PTM1 ,  PL ,  PI, &
 & PLUDE  , PLU   , PMFU   , PMFD ,&
 & PTENT  , PGTENT, PTENQ  , PGTENQ, &
 & PTENL  , PGTENL, PTENI  , PGTENI, &
 & PCLC   , PFPLSL , PFPLSN, &
 & PFHPSL , PFHPSN, PCOVPTOT )  

!**** *CLOUDSC2*  - COMPUTES CLOUD COVER, CLOUD LIQUID WATER/ICE
!                   AND LARGE-SCALE CONDENSATION.
!                   PROGNOSTIC CLOUD LIQUID WATER/ICE VARIABLES.
!                   (non-linear)

!     PURPOSE.
!     --------
!         THIS ROUTINE COMPUTES CLOUD AMOUNTS WHICH ARE REQUIRED BY THE
!     RADIATION SCHEME FOR COMPUTATION OF THE FLUXES AND HE PHYSICAL
!     TENDENCIES OF THE PROGNOSTIC VARIABLES T, Q, QL AND QI DUE TO
!     THE WATER PHASE CHANGES

!**   INTERFACE.
!     ----------
!              *CLOUDSC2* IS CALLED FROM *CALLPAR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                     S
!    *PQM1*         SPECIFIC HUMIDITY AT T-1                     KG/KG
!    *PQS*          SATURATION SPECIFIC HUMIDITY                 KG/KG
!    *PTM1*         TEMPERATURE AT T-1                             K
!    *PL*           CLOUD LIQUID WATER                           KG/KG       
!    *PI*           CLOUD ICE                                    KG/KG       
!    *PAPHP1*       PRESSURE AT HALF LEVELS AT T-1                PA    [#]
!    *PAPP1*        PROVISIONAL PRESSURE ON FULL LEVELS           PA
!    *PLUDE*        DETRAINED LIQUID WATER                       KG/(M3*S)
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
!    *PMFU*         MASS FLUX IN CONVECTIVE UPDRAFTS             KG/(M2*S)
!    *PMFD*         MASS FLUX IN CONVECTIVE DOWNDRAFTS           KG/(M2*S)
!    *PGTENT*       ACCUMULATED TEMPERATURE TENDENCY               K/S
!    *PGTENQ*       ACCUMULATED MOISTURE TENDENCY                KG/(KG S)
!    *PGTENL*       ACCUMULATED CLOUD LIQUID WATER TENDENCY      KG/(KG S)      
!    *PGTENI*       ACCUMULATED CLOUD ICE  TENDENCY              KG/(KG S)      

!    OUTPUT PARAMETERS (REAL):

!    *PTENT*        PROCESS TEMPERATURE TENDENCY                   K/S
!    *PTENQ*        PROCESS MOISTURE TENDENCY                    KG/(KG S)
!    *PTENL*        PROCESS CLOUD LIQUID WATER TENDENCY          KG/(KG S)      
!    *PTENI*        PROCESS CLOUD ICE  TENDENCY                  KG/(KG S)      
!    *PCLC*         CLOUD COVER OF INDIVIDUAL LAYER
!    *PFHPSL*       ENTHALPY FLUX DUE TO LARGE SCALE RAIN         J/(M2*S)
!    *PFHPSN*       ENTHALPY FLUX DUE TO LARGE SCALE SNOW         J/(M2*S)
!    *PFPLSL*       LARGE SCALE RAIN FLUX                        KG/(M2*S)
!    *PFPLSN*       LARGE SCALE SNOW FLUX                        KG/(M2*S)

!    *PCOVPTOT*     PRECIPITATION FRACTION IN EACH LAYER

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!     2009-04-08  P. Lopez

!     MODIFICATIONS.
!     --------------
!     F. Vana  2-Sep-2020: Introduced to CY46R1

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RETV, RG  ,RCPD     ,&
 & RLVTT    ,RLSTT    ,RLMLT    ,RTT     ,RD
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU   ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,RVTMP2 
USE YOECLD   , ONLY : TECLD
USE YOECLDP  , ONLY : TECLDP
USE YOEPHLI  , ONLY : TEPHLI
USE YOPHNC   , ONLY : TPHNC

IMPLICIT NONE

TYPE(TECLD)       ,INTENT(INOUT) :: YDECLD
TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TPHNC)       ,INTENT(INOUT) :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
LOGICAL           ,INTENT(IN)    :: LDRAIN1D 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHP1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPP1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLUDE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOVPTOT(KLON,KLEV)

!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------

INTEGER(KIND=JPIM) :: JK, JL

!=======================================================================
!     TUNABLE CONSTANTS (to be moved to include files later)
!=======================================================================

! zscal is a scale factor that linearly reduces the variance between 
! qv=qv-crit and qv-qsat
! 0 = No scaling
! 1 = full scaling (i.e. variance=0 when qv=qsat)

REAL(KIND=JPRB) :: ZSCAL=0.9_JPRB

!=======================================================================

REAL(KIND=JPRB) :: ZTP1(KLON,KLEV), ZQP1(KLON,KLEV), ZL(KLON,KLEV), ZI(KLON,KLEV)
REAL(KIND=JPRB) :: ZLUDE(KLON,KLEV)
REAL(KIND=JPRB) :: ZQC(KLON,KLEV), ZQLWC(KLON,KLEV), ZQIWC(KLON,KLEV)
REAL(KIND=JPRB) :: ZRFL(KLON), ZSFL(KLON), ZDP(KLON,KLEV)
REAL(KIND=JPRB) :: ZLSDCP(KLON,KLEV), ZLFDCP(KLON,KLEV), ZLVDCP(KLON,KLEV)
REAL(KIND=JPRB) :: ZRFLN(KLON), ZSFLN(KLON), ZGDP(KLON)
REAL(KIND=JPRB) :: ZDQDT(KLON), ZDTDT(KLON), ZDLDT(KLON), ZDIDT(KLON)
REAL(KIND=JPRB) :: ZQCRIT(KLON), ZCOVPCLR(KLON), ZCOVPTOT(KLON) 
REAL(KIND=JPRB) :: ZDQSDTEMP(KLON), ZCORQS(KLON), ZDTGDP(KLON)
REAL(KIND=JPRB) :: ZQOLD(KLON),ZPP(KLON),ZDQ(KLON), ZQLIM(KLON)
REAL(KIND=JPRB) :: ZSCALM(KLEV)
REAL(KIND=JPRB) :: ZRFREEZE(KLON,KLEV)
REAL(KIND=JPRB) :: ZCONDL(KLON,KLEV), ZCONDI(KLON,KLEV)
REAL(KIND=JPRB) :: ZEVAPR(KLON,KLEV), ZEVAPS(KLON,KLEV)
REAL(KIND=JPRB) :: ZQSAT(KLON), ZFOEEW(KLON), ZFWAT(KLON)

REAL(KIND=JPRB) :: ZOEALFA, ZCLDL, ZCLDI, ZLNEW, ZINEW, ZPRR, ZPRS, ZLCRIT, ZD
REAL(KIND=JPRB) :: ZCKCODTL, ZCKCODTI, ZCONS, ZSNMLT, ZZZ, ZRN, ZSN
REAL(KIND=JPRB) :: ZOEALFAW, ZFAC, ZFACW, ZFACI, ZESDP, ZCOR
REAL(KIND=JPRB) :: ZPRTOT, ZDPR, ZPRECLR, ZDR, ZDR2
REAL(KIND=JPRB) :: ZQE, ZBETA, ZB
REAL(KIND=JPRB) :: ZCONS2, ZCONS3, ZMELTP2, Z3ES, Z4ES, ZQMAX
REAL(KIND=JPRB) :: ZQPD, ZQCD, ZFWATR, ZRFREEZE2
REAL(KIND=JPRB) :: ZQT, ZQTMST

REAL(KIND=JPRB) :: ZRH1,ZRH2,ZRH3,ZETA2,ZDETA1,ZDETA2 
REAL(KIND=JPRB) :: ZTRPAUS(KLON)
REAL(KIND=JPRB) :: ZCRH2,ZETA3
REAL(KIND=JPRB) :: ZRHO,ZRODQSDP,ZLDCP,DTDZMO,ZDQSDZ,ZDQC 
REAL(KIND=JPRB) :: ZSUPSAT 
REAL(KIND=JPRB) :: ZFAC1,ZFAC2,ZFAC3,ZFAC4 

REAL(KIND=JPRB) :: ZEPS1,ZEPS2
INTEGER(KIND=JPIM) :: IK,ICALL

LOGICAL :: LLO1, LLO2, LLFLAG(KLON)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "cuadjtqs.intfb.h"

!     ------------------------------------------------------------------
#include "fcttre.func.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUDSC2',0,ZHOOK_HANDLE)
ASSOCIATE(CETA=>YDECLD%CETA, &
 & RCLCRIT=>YDECLDP%RCLCRIT, RKCONV=>YDECLDP%RKCONV, RLMIN=>YDECLDP%RLMIN, &
 & RPECONS=>YDECLDP%RPECONS, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, RLPTRC=>YDEPHLI%RLPTRC, &
 & LEVAPLS2=>YDPHNC%LEVAPLS2)

!*         1.     SET-UP INPUT QUANTITIES
!                 -----------------------

!*         1.1    Set-up tunning parameters

! set up constants required

ZCKCODTL=2.0_JPRB*RKCONV*PTSPHY
ZCKCODTI=5.0_JPRB*RKCONV*PTSPHY
ZCONS2 =1.0_JPRB/(PTSPHY*RG)
ZCONS3 =RLVTT/RCPD
ZMELTP2=RTT+2.0_JPRB
ZQTMST=1.0_JPRB/PTSPHY

ZQMAX=0.5_JPRB
ZEPS1=1.E-12_JPRB
ZEPS2=1.E-10_JPRB

!     --------------------------------------------------------------------

!*         2.1    COMPUTE CRITICAL RELATIVE HUMIDITY AND RELATIVE HUMIDITY
!                 --------------------------------------------------------

! first guess values for T, q, ql and qi

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP1(JL,JK)=PTM1(JL,JK) + PTSPHY*PGTENT(JL,JK)
    ZQP1(JL,JK)=PQM1(JL,JK) + PTSPHY*PGTENQ(JL,JK)
    ZL  (JL,JK)=PL  (JL,JK) + PTSPHY*PGTENL(JL,JK)
    ZI  (JL,JK)=PI  (JL,JK) + PTSPHY*PGTENI(JL,JK)
  ENDDO
ENDDO

DO JK=1,KLEV

! Parameter for cloud formation

  ZSCALM(JK)=ZSCAL*MAX((CETA(JK)-0.2_JPRB),ZEPS1)**(0.2_JPRB)

  DO JL=KIDIA,KFDIA

! thermodynamic constants

    ZDP(JL,JK)=PAPHP1(JL,JK+1)-PAPHP1(JL,JK)
    ZZZ=1.0_JPRB/(RCPD+RCPD*RVTMP2*ZQP1(JL,JK))
    ZLFDCP(JL,JK)=RLMLT*ZZZ
    ZLSDCP(JL,JK)=RLSTT*ZZZ
    ZLVDCP(JL,JK)=RLVTT*ZZZ
    LLFLAG(JL)=.TRUE.
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.2    INITIALIZATION OF CLOUD AND PRECIPITATION ARRAYS
!                 ------------------------------------------------

!       Clear cloud and freezing arrays

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PCLC(JL,JK)=0.0_JPRB
    ZQC(JL,JK)=0.0_JPRB
    ZQLWC(JL,JK)=0.0_JPRB
    ZQIWC(JL,JK)=0.0_JPRB
    ZRFREEZE(JL,JK)=0.0_JPRB
    ZCONDL(JL,JK)=0.0_JPRB
    ZCONDI(JL,JK)=0.0_JPRB
    ZEVAPR(JL,JK)=0.0_JPRB
    ZEVAPS(JL,JK)=0.0_JPRB
    PCOVPTOT(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!       Set to zero precipitation fluxes at the top

DO JL=KIDIA,KFDIA
  ZRFL(JL)=0.0_JPRB
  ZSFL(JL)=0.0_JPRB
  PFPLSL(JL,1)=0.0_JPRB
  PFPLSN(JL,1)=0.0_JPRB
  ZCOVPTOT(JL)=0.0_JPRB
  ZCOVPCLR(JL)=0.0_JPRB
ENDDO

! Eta value at tropopause
DO JL=KIDIA,KFDIA
  ZTRPAUS(JL)=0.1_JPRB
ENDDO
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    LLO1=CETA(JK) > 0.1_JPRB .AND. CETA(JK) < 0.4_JPRB .AND. &
       & ZTP1(JL,JK) > ZTP1(JL,JK+1)
    IF(LLO1) THEN
      ZTRPAUS(JL)=CETA(JK)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*        3. COMPUTE LAYER CLOUD AMOUNTS
!            ---------------------------

! Large loop over KLEV 
! Calculates 
!   1. diagnostic CC and QL
!   2. Convective CC and QL
!   3. Rainfall 

DO JK=KTDIA,KLEV

!       3.1   INITIALIZATION

  DO JL=KIDIA,KFDIA

!-----------------------------------
! calculate dqs/dT correction factor
!-----------------------------------

    IF (LPHYLIN .OR. LDRAIN1D) THEN
      ZOEALFAW=0.545_JPRB*(TANH(0.17_JPRB*(ZTP1(JL,JK)-RLPTRC))+1.0_JPRB)
      IF (ZTP1(JL,JK) < RTT) THEN
        ZFWAT(JL)=ZOEALFAW
        Z3ES=R3IES
        Z4ES=R4IES
      ELSE
        ZFWAT(JL)=1.0_JPRB
        Z3ES=R3LES
        Z4ES=R4LES
      ENDIF
      ZFOEEW(JL) = R2ES*EXP(Z3ES*(ZTP1(JL,JK)-RTT)/(ZTP1(JL,JK)-Z4ES))
      ZESDP = ZFOEEW(JL)/PAPP1(JL,JK)
      IF (ZESDP > ZQMAX) THEN
        ZESDP=ZQMAX
      ENDIF
    ELSE
      ZFWAT(JL)=FOEALFA(ZTP1(JL,JK))
      ZFOEEW(JL)=FOEEWM(ZTP1(JL,JK))
      ZESDP=ZFOEEW(JL)/PAPP1(JL,JK)
    ENDIF
    ZFACW=R5LES/((ZTP1(JL,JK)-R4LES)**2)
    ZFACI=R5IES/((ZTP1(JL,JK)-R4IES)**2)
    ZFAC=ZFWAT(JL)*ZFACW+(1.0_JPRB-ZFWAT(JL))*ZFACI
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
    ZDQSDTEMP(JL)=ZFAC*ZCOR*PQS(JL,JK)
    ZCORQS(JL)=1.0_JPRB+ZCONS3*ZDQSDTEMP(JL)

! use clipped state

    ZQLIM(JL)=ZQP1(JL,JK)
    IF (ZQP1(JL,JK)>PQS(JL,JK)) ZQLIM(JL)=PQS(JL,JK)

! set up critical value of humidity

    ZETA3=ZTRPAUS(JL)
    ZRH1=1.0_JPRB
    ZRH2=0.35_JPRB+0.14_JPRB*((ZETA3-0.25_JPRB)/0.15_JPRB)**2&
      & +0.04_JPRB*MIN(ZETA3-0.25_JPRB,0.0_JPRB)/0.15_JPRB
    ZRH3=1.0_JPRB
    ZDETA2=0.3_JPRB
    ZDETA1=0.09_JPRB+0.16_JPRB*(0.4_JPRB-ZETA3)/0.3_JPRB
    IF (CETA(JK) < ZETA3) THEN 
      ZCRH2=ZRH3
    ELSE IF (CETA(JK) >= ZETA3 .AND. CETA(JK) < (ZETA3+ZDETA2)) THEN 
      ZCRH2=ZRH3+(ZRH2-ZRH3)*((CETA(JK)-ZETA3)/ZDETA2)
    ELSE IF (CETA(JK) >= (ZETA3+ZDETA2) .AND. CETA(JK) < (1.0_JPRB-ZDETA1)) THEN 
      ZCRH2=ZRH2
    ELSE IF (CETA(JK) >= (1.0_JPRB-ZDETA1)) THEN 
      ZCRH2=ZRH1+(ZRH2-ZRH1)*((1.0_JPRB-CETA(JK))/ZDETA1)**0.5_JPRB
    ENDIF
! Allow ice supersaturation at cold temperatures
    IF (ZTP1(JL,JK) < RTICE) THEN
      ZSUPSAT=1.8_JPRB - 3.E-03_JPRB*ZTP1(JL,JK)
    ELSE
      ZSUPSAT=1.0_JPRB
    ENDIF
    ZQSAT(JL)=PQS(JL,JK)*ZSUPSAT
    ZQCRIT(JL)=ZCRH2*ZQSAT(JL)
  ENDDO

! Simple UNIFORM distribution of total water from Letreut & Li (90)

  DO JL=KIDIA,KFDIA
    ZQT=ZQP1(JL,JK)+ZL(JL,JK)+ZI(JL,JK)
    IF (ZQT <= ZQCRIT(JL)) THEN
      PCLC(JL,JK)=0.0_JPRB
      ZQC(JL,JK)=0.0_JPRB
    ELSEIF (ZQT >= ZQSAT(JL)) THEN
      PCLC(JL,JK)=1.0_JPRB
      ZQC(JL,JK)=(1.0_JPRB-ZSCALM(JK))*(ZQSAT(JL)-ZQCRIT(JL))
    ELSE
      ZQPD=ZQSAT(JL)-ZQT
      ZQCD=ZQSAT(JL)-ZQCRIT(JL)
      PCLC(JL,JK)=1.0_JPRB-SQRT(ZQPD/(ZQCD-ZSCALM(JK)*(ZQT-ZQCRIT(JL))))  
      ZQC(JL,JK)=(ZSCALM(JK)*ZQPD+(1.0_JPRB-ZSCALM(JK))*ZQCD) &
       & *PCLC(JL,JK)**2  
    ENDIF
  ENDDO

! Add convective component

  DO JL=KIDIA,KFDIA
    ZGDP(JL)=RG/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))
    ZLUDE(JL,JK)=PLUDE(JL,JK)*PTSPHY*ZGDP(JL)
    IF (JK<KLEV) THEN
      LLO1=ZLUDE(JL,JK)>=RLMIN.AND.PLU(JL,JK+1)>=ZEPS2
    ELSE
      LLO1=.FALSE.
    ENDIF
    IF (LLO1) THEN
      PCLC(JL,JK)=PCLC(JL,JK)+&
       & (1.0_JPRB-PCLC(JL,JK))*(1.0_JPRB-EXP(-ZLUDE(JL,JK)/PLU(JL,JK+1)))  
      ZQC(JL,JK)=ZQC(JL,JK)+ZLUDE(JL,JK)
    ENDIF
  ENDDO

! Add compensating subsidence component

  DO JL=KIDIA,KFDIA
    ZFAC1=1.0_JPRB/(RD*ZTP1(JL,JK))
    ZRHO = PAPP1(JL,JK)*ZFAC1
    ZFAC2=1.0_JPRB/(PAPP1(JL,JK)-RETV*ZFOEEW(JL))
    ZRODQSDP = -ZRHO*PQS(JL,JK)*ZFAC2
    ZLDCP = ZFWAT(JL)*ZLVDCP(JL,JK) + (1.0_JPRB-ZFWAT(JL))*ZLSDCP(JL,JK)
    ZFAC3=1.0_JPRB/(1.0_JPRB + ZLDCP*ZDQSDTEMP(JL))
    DTDZMO = RG*(1.0_JPRB/RCPD - ZLDCP*ZRODQSDP)*ZFAC3
    ZDQSDZ = ZDQSDTEMP(JL)*DTDZMO - RG*ZRODQSDP
    ZFAC4=1.0_JPRB/ZRHO
    ZDQC = MIN(ZDQSDZ*(PMFU(JL,JK)+PMFD(JL,JK))*PTSPHY*ZFAC4, ZQC(JL,JK))
    ZQC(JL,JK) = ZQC(JL,JK) - ZDQC
  ENDDO

! New cloud liquid/ice contents and condensation rates (liquid/ice)

  DO JL=KIDIA,KFDIA
    ZQLWC(JL,JK)=ZQC(JL,JK)*ZFWAT(JL)
    ZQIWC(JL,JK)=ZQC(JL,JK)*(1.0_JPRB-ZFWAT(JL))
    ZCONDL(JL,JK)= (ZQLWC(JL,JK)-ZL(JL,JK))*ZQTMST
    ZCONDI(JL,JK)= (ZQIWC(JL,JK)-ZI(JL,JK))*ZQTMST
  ENDDO


! Calculate precipitation overlap. 
! Simple form based on Maximum Overlap.

  DO JL=KIDIA,KFDIA
    IF (PCLC(JL,JK) > ZCOVPTOT(JL)) THEN ! total rain frac
      ZCOVPTOT(JL) = PCLC(JL,JK)
    ENDIF
    ZCOVPCLR(JL)=ZCOVPTOT(JL)-PCLC(JL,JK)   ! clear sky frac
    ZCOVPCLR(JL)=MAX(ZCOVPCLR(JL),0.0_JPRB)
  ENDDO

!*         3.3    CALCULATE PRECIPITATION

! Melting of incoming snow

  DO JL=KIDIA,KFDIA
    IF (ZSFL(JL) /= 0.0_JPRB) THEN
      ZCONS=ZCONS2*ZDP(JL,JK)/ZLFDCP(JL,JK)
      ZSNMLT=MIN(ZSFL(JL),ZCONS*MAX(0.0_JPRB,(ZTP1(JL,JK)-ZMELTP2)))
      ZRFLN(JL)=ZRFL(JL)+ZSNMLT
      ZSFLN(JL)=ZSFL(JL)-ZSNMLT
      ZTP1(JL,JK)=ZTP1(JL,JK)-ZSNMLT/ZCONS
    ELSE
      ZRFLN(JL)=ZRFL(JL)
      ZSFLN(JL)=ZSFL(JL)
    ENDIF
  ENDDO

  DO JL=KIDIA,KFDIA

!   Diagnostic calculation of rain production from cloud liquid water

    IF (PCLC(JL,JK)>ZEPS2) THEN
      IF (LEVAPLS2 .OR. LDRAIN1D) THEN
        ZLCRIT=1.9_JPRB*RCLCRIT
      ELSE
        ZLCRIT=RCLCRIT*2._JPRB
      ENDIF
      ZCLDL=ZQLWC(JL,JK)/PCLC(JL,JK) ! in-cloud liquid 
      ZD=ZCKCODTL*(1.0_JPRB-EXP(-(ZCLDL/ZLCRIT)**2))
      ZLNEW=PCLC(JL,JK)*ZCLDL*EXP(-ZD)
      ZPRR=ZQLWC(JL,JK)-ZLNEW
      ZQLWC(JL,JK)=ZQLWC(JL,JK)-ZPRR
    ELSE
      ZPRR=0.0_JPRB
    ENDIF

!   Diagnostic calculation of snow production from cloud ice

    IF (PCLC(JL,JK)>ZEPS2) THEN
      IF (LEVAPLS2 .OR. LDRAIN1D) THEN
        ZLCRIT=1.E-04_JPRB
      ELSE
        ZLCRIT=RCLCRIT*2._JPRB
      ENDIF
      ZCLDI=ZQIWC(JL,JK)/PCLC(JL,JK) ! in-cloud ice 
      ZD=ZCKCODTI*EXP(0.025_JPRB*(ZTP1(JL,JK)-RTT))*(1.0_JPRB-EXP(-(ZCLDI/ZLCRIT)**2))
      ZINEW=PCLC(JL,JK)*ZCLDI*EXP(-ZD)
      ZPRS=ZQIWC(JL,JK)-ZINEW
      ZQIWC(JL,JK)=ZQIWC(JL,JK)-ZPRS
    ELSE
      ZPRS=0.0_JPRB
    ENDIF

!   New precipitation (rain + snow)

    ZDR=ZCONS2*ZDP(JL,JK)*(ZPRR+ZPRS)

!   Rain fraction (different from cloud liquid water fraction!)

    IF (ZTP1(JL,JK) < RTT) THEN
      ZRFREEZE(JL,JK)=ZCONS2*ZDP(JL,JK)*ZPRR
      ZFWATR=0.0_JPRB
    ELSE
      ZFWATR=1.0_JPRB
    ENDIF 

    ZRN=ZFWATR*ZDR
    ZSN=(1.0_JPRB-ZFWATR)*ZDR
    ZRFLN(JL)=ZRFLN(JL)+ZRN
    ZSFLN(JL)=ZSFLN(JL)+ZSN

!   Precip evaporation

    ZPRTOT=ZRFLN(JL)+ZSFLN(JL)
    LLO2=ZPRTOT>ZEPS2 .AND. ZCOVPCLR(JL)>ZEPS2 .AND. (LEVAPLS2 .OR. LDRAIN1D)
    IF (LLO2) THEN 

      ZPRECLR=ZPRTOT*ZCOVPCLR(JL)/ZCOVPTOT(JL)

!     This is the humidity in the moistest zcovpclr region

      ZQE=PQS(JL,JK)-(PQS(JL,JK)-ZQLIM(JL))*ZCOVPCLR(JL)/&
       & (1.0_JPRB-PCLC(JL,JK))**2  
      ZBETA=RG*RPECONS*(SQRT(PAPP1(JL,JK) &
       & /PAPHP1(JL,KLEV+1))/5.09E-3_JPRB*ZPRECLR &
       & /ZCOVPCLR(JL))**0.5777_JPRB  

!     implicit solution:
      ZB=PTSPHY*ZBETA*(PQS(JL,JK)-ZQE)/(1.0_JPRB+ZBETA*PTSPHY*ZCORQS(JL))

!     exact solution:
!     ZB=(PQS(JL,JK)-ZQE)*(_ONE_-EXP(-ZBETA*ZCORQS(JL)*PTSPHY))/ZCORQS(JL)

      ZDTGDP(JL)=PTSPHY*RG/(PAPHP1(JL,JK+1)-PAPHP1(JL,JK))

      ZDPR=ZCOVPCLR(JL)*ZB/ZDTGDP(JL)
      ZDPR=MIN(ZDPR,ZPRECLR)
      ZPRECLR=ZPRECLR-ZDPR  ! take away from clr sky flux
      IF (ZPRECLR <= 0.0_JPRB) ZCOVPTOT(JL) = PCLC(JL,JK) !reset
      PCOVPTOT(JL,JK) = ZCOVPTOT(JL)

! warm proportion
      ZEVAPR(JL,JK)=ZDPR*ZRFLN(JL)/ZPRTOT
      ZRFLN(JL)=ZRFLN(JL)-ZEVAPR(JL,JK)

! ice proportion
      ZEVAPS(JL,JK)=ZDPR*ZSFLN(JL)/ZPRTOT
      ZSFLN(JL)=ZSFLN(JL)-ZEVAPS(JL,JK)
    ENDIF

  ENDDO

! Update of T and Q tendencies due to:
!  - condensation/evaporation of cloud liquid water/ice
!  - detrainment of convective cloud condensate
!  - evaporation of precipitation
!  - freezing of rain (impact on T only).

  DO JL=KIDIA,KFDIA
    ZDQDT(JL)=-(ZCONDL(JL,JK)+ZCONDI(JL,JK)) &
            & +(PLUDE(JL,JK)+ZEVAPR(JL,JK)+ZEVAPS(JL,JK))*ZGDP(JL)

    ZDTDT(JL)= ZLVDCP(JL,JK)*ZCONDL(JL,JK)+ZLSDCP(JL,JK)*ZCONDI(JL,JK) &
            & -( ZLVDCP(JL,JK)*ZEVAPR(JL,JK)+ZLSDCP(JL,JK)*ZEVAPS(JL,JK) &
            & +PLUDE(JL,JK)*(ZFWAT(JL)*ZLVDCP(JL,JK)+ &
            & (1.0_JPRB-ZFWAT(JL))*ZLSDCP(JL,JK))&
            & -(ZLSDCP(JL,JK)-ZLVDCP(JL,JK))*ZRFREEZE(JL,JK) )&
            & *ZGDP(JL)

! first guess T and Q
    ZTP1(JL,JK)=ZTP1(JL,JK) + PTSPHY*ZDTDT(JL)
    ZQP1(JL,JK)=ZQP1(JL,JK) + PTSPHY*ZDQDT(JL)

    ZPP(JL)=PAPP1(JL,JK)
    ZQOLD(JL)=ZQP1(JL,JK)
  ENDDO

! clipping of final qv

  IK=JK
  ICALL=0
  IF (LPHYLIN .OR. LDRAIN1D) THEN
    CALL CUADJTQS ( KIDIA, KFDIA, KLON, KLEV, IK,&
     & ZPP  , ZTP1  , ZQP1 , LLFLAG, ICALL  )  
  ELSE
    CALL CUADJTQ  ( YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK,&
     & ZPP  , ZTP1  , ZQP1 , LLFLAG, ICALL  )  
  ENDIF

  DO JL=KIDIA,KFDIA
    ZDQ(JL)=MAX(0.0_JPRB,ZQOLD(JL)-ZQP1(JL,JK))
    ZDR2=ZCONS2*ZDP(JL,JK)*ZDQ(JL)
! Update rain fraction and freezing.
! Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
    IF (ZTP1(JL,JK) < RTT) THEN
      ZRFREEZE2=ZFWAT(JL)*ZDR2
      ZFWATR=0.0_JPRB
    ELSE
      ZRFREEZE2=0.0_JPRB
      ZFWATR=1.0_JPRB
    ENDIF 
    ZRN=ZFWATR*ZDR2
    ZSN=(1.0_JPRB-ZFWATR)*ZDR2
! Note: The extra condensation due to the adjustment goes directly to precipitation
    ZCONDL(JL,JK)=ZCONDL(JL,JK)+ZFWATR*ZDQ(JL)*ZQTMST
    ZCONDI(JL,JK)=ZCONDI(JL,JK)+(1.0_JPRB-ZFWATR)*ZDQ(JL)*ZQTMST
    ZRFLN(JL)=ZRFLN(JL)+ZRN
    ZSFLN(JL)=ZSFLN(JL)+ZSN
    ZRFREEZE(JL,JK)=ZRFREEZE(JL,JK)+ZRFREEZE2
  ENDDO  

  DO JL=KIDIA,KFDIA
    ZDQDT(JL)=-(ZCONDL(JL,JK)+ZCONDI(JL,JK)) &
            & +(PLUDE(JL,JK)+ZEVAPR(JL,JK)+ZEVAPS(JL,JK))*ZGDP(JL)

    ZDTDT(JL)= ZLVDCP(JL,JK)*ZCONDL(JL,JK)+ZLSDCP(JL,JK)*ZCONDI(JL,JK) &
            & -( ZLVDCP(JL,JK)*ZEVAPR(JL,JK)+ZLSDCP(JL,JK)*ZEVAPS(JL,JK) &
            & +PLUDE(JL,JK)*(ZFWAT(JL)*ZLVDCP(JL,JK)+ &
            & (1.0_JPRB-ZFWAT(JL))*ZLSDCP(JL,JK))&
            & -(ZLSDCP(JL,JK)-ZLVDCP(JL,JK))*ZRFREEZE(JL,JK) )&
            & *ZGDP(JL)

    ZDLDT(JL)=(ZQLWC(JL,JK)-ZL(JL,JK))*ZQTMST

    ZDIDT(JL)=(ZQIWC(JL,JK)-ZI(JL,JK))*ZQTMST

    PTENQ(JL,JK)=ZDQDT(JL)
    PTENT(JL,JK)=ZDTDT(JL)
    PTENL(JL,JK)=ZDLDT(JL)
    PTENI(JL,JK)=ZDIDT(JL)

    PFPLSL(JL,JK+1)=ZRFLN(JL)
    PFPLSN(JL,JK+1)=ZSFLN(JL)
  ENDDO

! record rain flux for next level

  DO JL=KIDIA,KFDIA
    ZRFL(JL)=ZRFLN(JL)
    ZSFL(JL)=ZSFLN(JL)
  ENDDO

ENDDO  !jk

!*     ENTHALPY FLUXES DUE TO PRECIPITATION
!      ------------------------------------

DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    PFHPSL(JL,JK)=-PFPLSL(JL,JK)*RLVTT
    PFHPSN(JL,JK)=-PFPLSN(JL,JK)*RLSTT
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUDSC2',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUDSC2

