! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLOUDST ( YDEPHLI,YDPHNC,YDECLD,YDECLDP,KIDIA , KFDIA , KLON   , KTDIA , KLEV   , LDRAIN1D,&
 & PTSPHY,&
 & PAPHP1 , PAPP1 ,  PQM1  , PQS   , PTM1,&
 & PLUDE  , PLU,&
 & PTENT  , PGTENT , PTENQ ,PGTENQ , &
 & PCLC   , PQIWP  , PQLWP ,PFPLSL , PFPLSN,&
 & PFHPSL , PFHPSN,  PCOVPTOT )  

!**** *CLOUDST*  - COMPUTES CLOUD COVER AND LIQUID WATER PATH
!                  AND LARGE-SCALE CONDENSATION

!     PURPOSE.
!     --------
!         THIS ROUTINE COMPUTES CLOUD AMOUNTS WHICH ARE REQUIRED BY THE
!     RADIATION SCHEME FOR COMPUTATION OF THE FLUXES AND HE PHYSICAL
!     TENDENCIES OF THE TWO PROGNOSTIC VARIABLES T AND Q DUE TO
!     THE WATER PHASE CHANGES

!**   INTERFACE.
!     ----------
!              *CLOUDST* IS CALLED FROM *CALLPAR*

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
!    *PAPHP1*       PRESSURE AT HALF LEVELS AT T-1                PA    [#]
!    *PAPP1*        PROVISIONAL PRESSURE ON FULL LEVELS           PA
!    *PLUDE*        DETRAINED LIQUID WATER                       KG/(M3*S)
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG

!    UPDATED PARAMETERS (REAL):
!    *PTENT*        TEMPERATURE TENDENCY                          K/S
!    *PGTENT*       TEMPERATURE TENDENCY TO DEFINE UPDATED STATE  K/S
!    *PTENQ*        MOISTURE TENDENCY                            KG/(KG S)
!    *PGTENQ*       MOISTURE TENDENCY TO DEFINE UPDATED STATE    KG/(KG S)

!    OUTPUT PARAMETERS (REAL):

!    *PCLC*         CLOUD COVER OF INDIVIDUAL LAYER
!    *PQIWP*        ICE WATER CONTENT                             KG/KG
!    *PQLWP*        LIQUID WATER CONTENT                          KG/KG

!    *PFHPSL*       ENTHALPY FLUX DUE TO LARGE SCALE RAIN         J/(M2*S)
!    *PFHPSN*       ENTHALPY FLUX DUE TO LARGE SCALE SNOW         J/(M2*S)
!    *PFPLSL*       LARGE SCALE RAIN FLUX                        KG/(M2*S)
!    *PFPLSN*       LARGE SCALE SNOW FLUX                        KG/(M2*S)

!    *PCOVPTOT*     PRECIPITATION FRACTION IN EACH LAYER

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!     2002-09-12  A. Tompkins and M. Janiskova

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Lopez      11-01-2007  Modified usage of LDRAIN
!        P.Lopez       22-Oct-2007 Retuning of autoconversion coefficient
!        A. Geer       01-Oct-2008 Name change to LDRAIN1D to reflect usage
!        A. Geer       13-Apr-2010 Output 2D precip fraction for all-sky radiances
!        M. Janiskova  11-Mar-2011 Cleaning and simplification of AG modification
!        F. Vana       18-May-2012 Cleaning
!        F.Vana        11-Sep-2020 Consistent thermodynamics
!        P.Lopez       15-Jan-2021 Revised critical relative humidity and reduced evaporation

!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RETV, RG  ,RCPD     ,&
 &                    RLVTT    ,RLSTT    ,RLMLT    ,RTT  
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU   ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,RVTMP2 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLUDE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQIWP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQLWP(KLON,KLEV) 
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

REAL(KIND=JPRB) :: ZTP1(KLON,KLEV), ZQP1(KLON,KLEV)
REAL(KIND=JPRB) :: ZLUDE(KLON,KLEV), ZQLWC(KLON,KLEV)
REAL(KIND=JPRB) :: ZRFL(KLON), ZSFL(KLON), ZDP(KLON,KLEV), ZDR(KLON)
REAL(KIND=JPRB) :: ZLSDCP(KLON,KLEV), ZLFDCP(KLON,KLEV), ZLVDCP(KLON,KLEV)
REAL(KIND=JPRB) :: ZRFLN(KLON), ZSFLN(KLON),ZGDP(KLON)
REAL(KIND=JPRB) :: ZFWAT(KLON,KLEV), ZDQDT(KLON), ZDTDT(KLON)
REAL(KIND=JPRB) :: ZQCRIT(KLON), ZCOVPTOT(KLON), ZCOVPCLR(KLON)
REAL(KIND=JPRB) :: ZDQSDTEMP(KLON), ZCORQS(KLON), ZDTGDP(KLON)
REAL(KIND=JPRB) :: ZQOLD(KLON),ZPP(KLON),ZDQ(KLON), ZQLIM(KLON)
REAL(KIND=JPRB) :: ZRFREEZE(KLON,KLEV)

REAL(KIND=JPRB) :: ZCRH(KLON,KLEV),ZSCALM(KLEV)
REAL(KIND=JPRB) :: ZTRPAUS(KLON)

REAL(KIND=JPRB) :: ZOEALFA, ZCLD, ZLNEW, ZPR, ZLCRIT, ZD
REAL(KIND=JPRB) :: ZCKCODT, ZCONS, ZSNMLT, ZZZ, ZRN, ZSN
REAL(KIND=JPRB) :: ZOEALFAW, ZALFAW, ZFAC, ZFACW, ZFACI, ZESDP, ZCOR
REAL(KIND=JPRB) :: ZFOEEW, ZPRTOT, ZDPR, ZDPRSFL, ZDPSSFL, ZPRECLR
REAL(KIND=JPRB) :: ZQE, ZBETA, ZB
REAL(KIND=JPRB) :: ZCONS2, ZCONS3, ZMELTP2, Z3ES, Z4ES, ZQMAX
REAL(KIND=JPRB) :: ZQPD, ZQCD, ZFWATR, ZRFREEZE2

REAL(KIND=JPRB) :: ZEPS1, ZEPS2
INTEGER(KIND=JPIM) :: IK, ICALL
REAL(KIND=JPRB) :: ZETA1, ZETA3, ZRH1, ZRH2, ZRH3, ZDETA1, ZPECONS

LOGICAL :: LLO1, LLO2, LLO3, LLFLAG(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "cuadjtqs.intfb.h"

!     ------------------------------------------------------------------
#include "fcttre.func.h"
!     ------------------------------------------------------------------

!*         1.     SET-UP INPUT QUANTITIES
!                 -----------------------

!*         1.1    Set-up tunning parameters

IF (LHOOK) CALL DR_HOOK('CLOUDST',0,ZHOOK_HANDLE)
ASSOCIATE(CETA=>YDECLD%CETA, &
 & RCLCRIT=>YDECLDP%RCLCRIT, RKCONV=>YDECLDP%RKCONV, RLMIN=>YDECLDP%RLMIN, &
 & RPECONS=>YDECLDP%RPECONS, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, RLPTRC=>YDEPHLI%RLPTRC, &
 & LEVAPLS2=>YDPHNC%LEVAPLS2)
 
! set up constants required

ZCKCODT=2.0_JPRB*RKCONV*PTSPHY
ZCONS2 =1.0_JPRB/(PTSPHY*RG)
ZCONS3 =2.5E6_JPRB/1005._JPRB
ZMELTP2=RTT+2.0_JPRB
ZPECONS=1.E-5_JPRB

ZQMAX=0.5_JPRB
ZEPS1=1.E-12_JPRB
ZEPS2=1.E-10_JPRB

! first guess values for T and q
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP1(JL,JK)=PTM1(JL,JK) + PTSPHY*PGTENT(JL,JK)
    ZQP1(JL,JK)=PQM1(JL,JK) + PTSPHY*PGTENQ(JL,JK)
  ENDDO
ENDDO

! Eta value at tropopause.
DO JL=KIDIA,KFDIA
  ZTRPAUS(JL)=0.1_JPRB
ENDDO
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    LLO3=CETA(JK) > 0.1_JPRB .AND. CETA(JK) < 0.4_JPRB .AND. &
       & ZTP1(JL,JK) > ZTP1(JL,JK+1)
    IF(LLO3) THEN
      ZTRPAUS(JL)=CETA(JK)
    ENDIF
  ENDDO
ENDDO

!     --------------------------------------------------------------------

!*         2.1    COMPUTE CRITICAL RELATIVE HUMIDITY AND RELATIVE HUMIDITY
!                 --------------------------------------------------------

! first guess values for T and q
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP1(JL,JK)=PTM1(JL,JK) + PTSPHY*PGTENT(JL,JK)
    ZQP1(JL,JK)=PQM1(JL,JK) + PTSPHY*PGTENQ(JL,JK)
  ENDDO
ENDDO

DO JK=1,KLEV

! Critical relative humidity

  DO JL=KIDIA,KFDIA
    ZETA3 = ZTRPAUS(JL)
    ZRH1 = 0.9_JPRB
    ZRH2 = 0.5_JPRB
    ZRH3 = 1.0_JPRB
    ZETA1 = 0.75_JPRB
    ZDETA1 = 1.0_JPRB - ZETA1
    IF (CETA(JK) < ZETA3) THEN 
      ZCRH(JL,JK) = ZRH3 + (ZRH2 - ZRH3) * CETA(JK) / ZETA3
    ELSEIF (CETA(JK) >= ZETA3 .AND. CETA(JK) < ZETA1) THEN 
      ZCRH(JL,JK) = ZRH2
    ELSEIF (CETA(JK) >= ZETA1) THEN 
      ZCRH(JL,JK) = ZRH1 + (ZRH2 - ZRH1) * ((1.0_JPRB - CETA(JK)) / ZDETA1)**0.5_JPRB
    ENDIF
  ENDDO 

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
    PQIWP(JL,JK)=0.0_JPRB
    PQLWP(JL,JK)=0.0_JPRB
    ZQLWC(JL,JK)=0.0_JPRB
    ZRFREEZE(JL,JK)=0.0_JPRB
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

    ZDR(JL)=0.0_JPRB

!-----------------------------------
! calculate dqs/dT correction factor
!-----------------------------------

    IF (LPHYLIN .OR. LDRAIN1D) THEN
      ZOEALFAW=0.545_JPRB*(TANH(0.17_JPRB*(ZTP1(JL,JK)-RLPTRC))+1.0_JPRB)
      IF (ZTP1(JL,JK) < RTT) THEN
        ZALFAW=ZOEALFAW
        Z3ES=R3IES
        Z4ES=R4IES
      ELSE
        ZALFAW=1.0_JPRB
        Z3ES=R3LES
        Z4ES=R4LES
      ENDIF
      ZFOEEW = R2ES*EXP(Z3ES*(ZTP1(JL,JK)-RTT)/(ZTP1(JL,JK)-Z4ES))
      ZESDP = ZFOEEW/PAPP1(JL,JK)
      IF (ZESDP > ZQMAX) THEN
        ZESDP=ZQMAX
      ENDIF
    ELSE
      ZALFAW=FOEALFA(ZTP1(JL,JK))
      ZESDP=FOEEWM(ZTP1(JL,JK))/PAPP1(JL,JK)
    ENDIF
    ZFACW=R5LES/((ZTP1(JL,JK)-R4LES)**2)
    ZFACI=R5IES/((ZTP1(JL,JK)-R4IES)**2)
    ZFAC=ZALFAW*ZFACW+(1.0_JPRB-ZALFAW)*ZFACI
    ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
    ZDQSDTEMP(JL)=ZFAC*ZCOR*PQS(JL,JK)
    ZCORQS(JL)=1.0_JPRB+ZCONS3*ZDQSDTEMP(JL)

! use clipped state

    ZQLIM(JL)=ZQP1(JL,JK)
    IF (ZQP1(JL,JK)>PQS(JL,JK)) ZQLIM(JL)=PQS(JL,JK)

! set up critical value of humidity

    ZQCRIT(JL)=ZCRH(JL,JK)*PQS(JL,JK)
  ENDDO

! Replace BETA distribution with simple UNIFORM from Letreut & Li (90)

  DO JL=KIDIA,KFDIA
    IF (ZQP1(JL,JK)<=ZQCRIT(JL)) THEN
      PCLC(JL,JK)=0.0_JPRB
      ZQLWC(JL,JK)=0.0_JPRB
    ELSEIF (ZQP1(JL,JK)>=PQS(JL,JK)) THEN
      PCLC(JL,JK)=1.0_JPRB
      ZQLWC(JL,JK)=(1.0_JPRB-ZSCALM(JK))*(PQS(JL,JK)-ZQCRIT(JL))
    ELSE
      ZQPD=PQS(JL,JK)-ZQP1(JL,JK)
      ZQCD=PQS(JL,JK)-ZQCRIT(JL)
      PCLC(JL,JK)=1.0_JPRB-SQRT(ZQPD/(ZQCD-ZSCALM(JK)*(ZQP1(JL,JK) &
       & -ZQCRIT(JL))))  
      ZQLWC(JL,JK)=(ZSCALM(JK)*ZQPD+(1.0_JPRB-ZSCALM(JK))*ZQCD) &
       & *PCLC(JL,JK)**2  
    ENDIF 
  ENDDO

! Add convective component

!DEC$ IVDEP
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
      ZQLWC(JL,JK)=ZQLWC(JL,JK)+ZLUDE(JL,JK)
    ENDIF
  ENDDO

! Calculate precipitation overlap. 
! Simple form based on Maximum Overlap.

  DO JL=KIDIA,KFDIA
    ZCOVPTOT(JL)=MAX(ZCOVPTOT(JL),PCLC(JL,JK)) ! total rain frac
    ZCOVPCLR(JL)=ZCOVPTOT(JL)-PCLC(JL,JK)      ! clear sky frac
    ZCOVPCLR(JL)=MAX(ZCOVPCLR(JL),0.0_JPRB)
  ENDDO

!*         3.3    CALCULATE PRECIPITATION

!    Melting of incoming snow

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

!    Diagnostic calculation of rain production

!DEC$ IVDEP
  DO JL=KIDIA,KFDIA
    IF (PCLC(JL,JK)>ZEPS2) THEN
      IF (LEVAPLS2 .OR. LDRAIN1D) THEN
        ZLCRIT=RCLCRIT
      ELSE
        ZLCRIT=RCLCRIT*2._JPRB
      ENDIF
      ZCLD=ZQLWC(JL,JK)/PCLC(JL,JK) ! in-cloud liquid 
      ZD=ZCKCODT*(1.0_JPRB-EXP(-(ZCLD/ZLCRIT)**2))
      ZLNEW=PCLC(JL,JK)*ZCLD*EXP(-ZD)
      ZPR=ZQLWC(JL,JK)-ZLNEW
      ZQLWC(JL,JK)=ZQLWC(JL,JK)-ZPR
    ELSE
      ZPR=0.0_JPRB
    ENDIF

!    New precipitation

    ZDR(JL)=ZDR(JL)+ZCONS2*ZDP(JL,JK)*ZPR
    IF (LPHYLIN .OR. LDRAIN1D) THEN
      ZOEALFA=0.545_JPRB*(TANH(0.17_JPRB*(ZTP1(JL,JK)-RLPTRC))+1.0_JPRB)
      IF (ZTP1(JL,JK) < RTT) THEN
        ZFWAT(JL,JK)=ZOEALFA
      ELSE
        ZFWAT(JL,JK)=1.0_JPRB
      ENDIF
    ELSE
      ZFWAT(JL,JK)=FOEALFA(ZTP1(JL,JK))
    ENDIF

! Rain fraction (different from cloud liquid water fraction!)
    IF (ZTP1(JL,JK) < RTT) THEN
      ZRFREEZE(JL,JK)=ZFWAT(JL,JK)*ZDR(JL)
      ZFWATR=0.0_JPRB
    ELSE
      ZFWATR=1.0_JPRB
    ENDIF 

    ZRN=ZFWATR*ZDR(JL)
    ZSN=(1.0_JPRB-ZFWATR)*ZDR(JL)
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
      ZBETA=RG*ZPECONS*(SQRT(PAPP1(JL,JK) &
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
      IF (ZPRECLR <= 0.0_JPRB) ZCOVPTOT(JL)=PCLC(JL,JK) !reset
      PCOVPTOT(JL,JK) = ZCOVPTOT(JL)

! warm proportion
      ZDPRSFL=ZDPR*ZRFLN(JL)/ZPRTOT
      ZRFLN(JL)=ZRFLN(JL)-ZDPRSFL

! ice proportion
      ZDPSSFL=ZDPR*ZSFLN(JL)/ZPRTOT
      ZSFLN(JL)=ZSFLN(JL)-ZDPSSFL
    ENDIF

!   Partitioning for cloud liquid and ice water content

    PQLWP(JL,JK)=ZFWAT(JL,JK)*ZQLWC(JL,JK)
    PQIWP(JL,JK)=(1.0_JPRB-ZFWAT(JL,JK))*ZQLWC(JL,JK)
  ENDDO

!   Incrementation of T and Q and fluxes' swap

  DO JL=KIDIA,KFDIA
    ZDQDT(JL)=-((ZRFLN(JL)-ZRFL(JL))+(ZSFLN(JL)-ZSFL(JL)))&
     & *(RG/ZDP(JL,JK))+PLUDE(JL,JK)*ZGDP(JL)  
    ZDTDT(JL)=(ZLVDCP(JL,JK)*(ZRFLN(JL)-ZRFL(JL))+ZLSDCP(JL,JK)&
     & *(ZSFLN(JL)-ZSFL(JL)))*(RG/ZDP(JL,JK))&
     & -PLUDE(JL,JK)*ZGDP(JL)* &
     & (ZFWAT(JL,JK)*ZLVDCP(JL,JK)+ &
     & (1.0_JPRB-ZFWAT(JL,JK))*ZLSDCP(JL,JK))&
     & +(ZLSDCP(JL,JK)-ZLVDCP(JL,JK))*ZRFREEZE(JL,JK)*ZGDP(JL)

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
    CALL CUADJTQ(YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK,&
     & ZPP, ZTP1, ZQP1, LLFLAG, ICALL)
  ENDIF

  DO JL=KIDIA,KFDIA
    ZDQ(JL)=MAX(0.0_JPRB,ZQOLD(JL)-ZQP1(JL,JK))
    ZDR(JL)=ZCONS2*ZDP(JL,JK)*ZDQ(JL) 
! Update rain fraction and freezing.
! Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
    IF (ZTP1(JL,JK) < RTT) THEN
      ZRFREEZE2=ZFWAT(JL,JK)*ZDR(JL)
      ZFWATR=0.0_JPRB
    ELSE
      ZRFREEZE2=0.0_JPRB
      ZFWATR=1.0_JPRB
    ENDIF 
    ZRN=ZFWATR*ZDR(JL)
    ZSN=(1.0_JPRB-ZFWATR)*ZDR(JL)
    ZRFLN(JL)=ZRFLN(JL)+ZRN
    ZSFLN(JL)=ZSFLN(JL)+ZSN
    ZRFREEZE(JL,JK)=ZRFREEZE(JL,JK)+ZRFREEZE2
  ENDDO  

  DO JL=KIDIA,KFDIA
    ZDQDT(JL)=-((ZRFLN(JL)-ZRFL(JL))+(ZSFLN(JL)-ZSFL(JL)))&
     & *(RG/ZDP(JL,JK))+PLUDE(JL,JK)*ZGDP(JL)  
    ZDTDT(JL)=(ZLVDCP(JL,JK)*(ZRFLN(JL)-ZRFL(JL))+ZLSDCP(JL,JK)&
     & *(ZSFLN(JL)-ZSFL(JL)))*(RG/ZDP(JL,JK))&
     & -PLUDE(JL,JK)*ZGDP(JL)* &
     & (ZFWAT(JL,JK)*ZLVDCP(JL,JK)+ &
     & (1.0_JPRB-ZFWAT(JL,JK))*ZLSDCP(JL,JK))& 
     & +(ZLSDCP(JL,JK)-ZLVDCP(JL,JK))*ZRFREEZE(JL,JK)*ZGDP(JL)

    PTENQ(JL,JK)=ZDQDT(JL)
    PTENT(JL,JK)=ZDTDT(JL)

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
IF (LHOOK) CALL DR_HOOK('CLOUDST',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUDST

