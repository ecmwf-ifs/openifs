! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUCALLN2 &
 & (  YDERAD,YDML_PHY_SLIN,YDML_PHY_EC,&
 & KIDIA,    KFDIA,    KLON,   KLEV,&
 & LDLAND, LDSLPHY, LDRAIN1D, &
 & PTSPHY,PVDIFTS,&
 & PTM1,     PQM1,     PUM1,     PVM1,&
 & PVERVEL,  PQHFL,    PAHFS,    PAPHM1,&
 & PAP,      PAPH,     PGEO,     PGEOH, PGAW,&
 & PTENT,    PTTENT,   PTENQ,    PTTENQ, &
 & PTENU,    PTTENU,   PTENV,    PTTENV, PARPRC,&
 & KTOPC,    KBASEC,   KTYPE,&
 & KCBOT,    KCTOP,    KBOTSC,   LDCUM,   LDSC,&
 & PLU,      PLUDE,    PMFU,     PMFD,&
 & PDIFCQ,   PDIFCS,   PFHPCL,   PFHPCN,&
 & PFPLCL,   PFPLCN,   PSTRCU,   PSTRCV,  PFCQLF ,PFCQIF,&
 & PMFUDE_RATE ,       PMFDDE_RATE ,      PCAPE,&
 & KTRAC,    PCM1,     PTENC,    PSCAV )  

!          *CUCALLN2* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                      *CUMASTRN2* (CUMULUS PARAMETERIZATION)
!                      *CUCCDIA* (CUMULUS CLOUDS FOR RADIATION)
!                      *CUSTRAT* (PBL_STRATOCUMULUS)

!           ***  SIMPLIFIED CONVECTION SCHEME ***

!           P. LOPEZ     E.C.M.W.F.     15/12/2003
!           Duplicated from M. TIEDTKE (1989)

!**   PURPOSE.
!     --------

!          *CUCALLN2* - INTERFACE FOR *CUMASTRN2*,*CUCCDIA* AND *CUSTRAT*:
!                       PROVIDES INPUT FOR CUMASTRN2, CUCCDIA AND CUSTRAT.
!                       RECEIVES UPDATED TENDENCIES, PRECIPITATION
!                       AND CONVECTIVE CLOUD PARAMETERS FOR RADIATION.

!**   INTERFACE.
!     ----------

!          *CUCALLN2* IS CALLED FROM *CALLPAR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                     S
!    *PTM1*         TEMPERATURE (T-1)                             K
!    *PQM1*         SPECIFIC HUMIDITY (T-1)                       KG/KG
!    *PUM1*         X-VELOCITY COMPONENT (T-1)                    M/S
!    *PVM1*         Y-VELOCITY COMPONENT (T-1)                    M/S
!    *PCM1*         CHEMICAL TRACERS (T-1)                        KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)     KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PAPHM1*       PRESSURE ON HALF LEVELS                       PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS           PA 
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PGAW*       NORMALISED GAUSSIAN QUADRATURE WEIGHT / NUMBER OF LONGITUDE POINTS
!                           LOCAL SUB-AREA == 4*RPI*RA**2 * PGAW
!    *PSCAV*        SCAVENGING COEFFICIENT                     UNITLESS


!    UPDATED PARAMETERS (REAL):

!   output tendencies relevat to the convection
!    *PTENT*        TEMPERATURE TENDENCY                         K/S 
!    *PTENQ*        MOISTURE TENDENCY                            KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                  M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                  M/S2 

!   input tendencies:
!    i/ The one cumulating tendencies for the output total tendecy 
!    *PTTENT*       TEMPERATURE TENDENCY                         K/S 
!    *PTTENQ*       MOISTURE TENDENCY                            KG/(KG S)
!    *PTTENU*       TENDENCY OF U-COMP. OF WIND (guess)          M/S2
!    *PTTENV*       TENDENCY OF V-COMP. OF WIND (guess)          M/S2 

!    *PTENC*        TENDENCY OF CHEMICAL TRACERS                 1/S
!    *PARPRC*       ACCUMULATED PRECIPITATION AMMOUNT           KG/(M2*S)
!                   FOR RADIATION CALCULATION

!    UPDATED PARAMETERS (INTEGER):

!    *KTOPC*        UPDATED CONVECTIVE CLOUD TOP LEVEL FOR RADIATION
!    *KBASEC*       UPDATED CONVECTIVE CLOUD BASE LEVEL FOR RADIATION

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION 
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (REAL):

!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
!    *PLUDE*        DETRAINED LIQUID WATER                      KG/(M2*S)
!    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                         KG/(M2*S)
!    *PDIFCQ*       CONVECTIVE FLUX OF SPECIFIC HUMIDITY        KG/(M2*S)
!    *PDIFCS*       CONVECTIVE FLUX OF HEAT                      J/(M2*S)
!    *PFHPCL*       ENTHALPY FLUX DUE TO CONVECTIVE RAIN         J/(M2*S)
!    *PFHPCN*       ENTHALPY FLUX DUE TO CONVECTIVE SNOW         J/(M2*S)
!    *PFPLCL*       CONVECTIVE RAIN FLUX                        KG/(M2*S)
!    *PFPLCN*       CONVECTIVE SNOW FLUX                        KG/(M2*S)
!    *PSTRTU*       CONVECTIVE FLUX OF U-MOMEMTUM         (M/S)*KG/(M2*S)
!    *PSTRTV*       CONVECTIVE FLUX OF V-MOMEMTUM         (M/S)*KG/(M2*S)
!    *PFCQLF*       CONVECTIVE CONDENSATION FLUX LIQUID         KG/(M2*S)
!    *PFCQIF*       CONVECTIVE CONDENSATION FLUX ICE            KG/(M2*S)

!               Diagnostics:
!    *PMFUDE_RATE* UPDRAFT DETRAINMENT RATE                     KG/(M3*S)
!    *PMFDDE_RATE* DOWNDRAFT DETRAINMENT RATE                   KG/(M3*S)
!    *PCAPE*       MAXIMUM pseudoadiabat CAPE                   (J/KG) 

!     EXTERNALS.
!     ----------

!          CUMASTRN2
!          CUCCDIA
!          CUSTRAT
!          SATUR
!          SATUR_1D

!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Lopez      14-02-2006  Added transport of tracers by mass flux
!        P. Lopez      11-01-2007  Added routine SATUR_1D for 1D-Var rain
!        A. Geer       10-01-2008  LDRAIN1D name change to reflect usage
!        F. Vana       18-May-2012 Cleaning
!        F. Vana       Oct-2013    Bug fix maintaining symetry with the NL scheme
!        P. Lopez      10-05-2016  Changed tendency computations to match full scheme

!-------------------------------------------------------------------------

USE MODEL_PHYSICS_ECMWF_MOD      , ONLY : MODEL_PHYSICS_ECMWF_TYPE
USE MODEL_PHYSICS_SIMPLINEAR_MOD , ONLY : MODEL_PHYSICS_SIMPLINEAR_TYPE
USE YOERAD                       , ONLY : TERAD
USE PARKIND1                     , ONLY : JPIM     ,JPRB
USE YOMHOOK                      , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST                       , ONLY : RG       ,RV       ,RCPD     ,RCPV     ,&
 &                                        RETV     ,RCW      ,RCS      ,RLVTT    ,RLSTT    ,&
 &                                        RTT      ,RALPW    ,RBETW    ,RGAMW    ,RALPS    ,&
 &                                        RBETS    ,RGAMS    ,RALPD    ,RBETD    ,RGAMD  
USE YOETHF                       , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                                        R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
 &                                        R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
 &                                        RTICECU, RTWAT_RTICECU_R, RTWAT_RTICE_R  
USE YOMCT3                       , ONLY : NSTEP

IMPLICIT NONE

TYPE(TERAD)                        ,INTENT(INOUT) :: YDERAD
TYPE(MODEL_PHYSICS_ECMWF_TYPE)     ,INTENT(INOUT) :: YDML_PHY_EC
TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE),INTENT(INOUT) :: YDML_PHY_SLIN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON)
LOGICAL           ,INTENT(IN)    :: LDSLPHY
LOGICAL           ,INTENT(IN)    :: LDRAIN1D 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTTENT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTTENQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTTENU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTTENV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCAV(KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PARPRC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTOPC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KBASEC(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KTYPE(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCBOT(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCTOP(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON)
LOGICAL           ,INTENT(OUT)   :: LDCUM(KLON)
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLUDE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCQ(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCL(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCN(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCV(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQLF(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQIF(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFUDE_RATE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFDDE_RATE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON)

!       3D DIAGNOSTICS FOR ERA40
!       OTHER CONVECTION DIAGNOSTICS
REAL(KIND=JPRB) ::     ZTP1(KLON,KLEV),        ZQP1(KLON,KLEV),&
 & ZUP1(KLON,KLEV),        ZVP1(KLON,KLEV),&
 & ZTU(KLON,KLEV),         ZQU(KLON,KLEV),&
 & ZQSAT(KLON,KLEV),       ZRAIN(KLON)         

REAL(KIND=JPRB) :: ZCP1(KLON,KLEV,KTRAC)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZENTHD(KLON,KLEV),ZENTHS(KLON,KLEV)  
REAL(KIND=JPRB) :: ZCONDFLL(KLON),ZCONDFLN(KLON)

INTEGER(KIND=JPIM) :: IFLAG, JK, JL, JN

REAL(KIND=JPRB) :: ZALFAW, ZCP, ZGDPH

LOGICAL :: LLTEST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuccdia.intfb.h"
#include "cumastrn2.intfb.h"
#include "custrat.intfb.h"
#include "satur.intfb.h"
#include "satur_1d.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "fcttrm.func.h"

!-----------------------------------------------------------------------

!*    0.1          STORE T,Q,X,Y TENDENCIES FOR FLUX COMPUTATIONS
!*                 ----------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCALLN2',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT22=>YDML_PHY_SLIN%YRECUMF2%NJKT22, &
 & LPHYLIN=>YDML_PHY_SLIN%YREPHLI%LPHYLIN, RLPTRC=>YDML_PHY_SLIN%YREPHLI%RLPTRC, &
 & LEPCLD=>YDML_PHY_EC%YREPHY%LEPCLD, LMFTRAC=>YDML_PHY_EC%YREPHY%LMFTRAC, &
 & LENCLD2=>YDML_PHY_SLIN%YRPHNC%LENCLD2,LEPCLD2=>YDML_PHY_SLIN%YRPHNC%LEPCLD2)
! Setup of tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENQ(JL,JK)=PTTENQ(JL,JK)
    PTENT(JL,JK)=PTTENT(JL,JK)
    PTENV(JL,JK)=PTTENV(JL,JK)
    PTENU(JL,JK)=PTTENU(JL,JK)
    ZENTHD(JL,JK)=0.0_JPRB
    ZENTHS(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    PFPLCL(JL,JK)=0.0_JPRB
    PFPLCN(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZUP1(JL,JK)=PUM1(JL,JK)+PTTENU(JL,JK)*PTSPHY
    ZVP1(JL,JK)=PVM1(JL,JK)+PTTENV(JL,JK)*PTSPHY
    ZTP1(JL,JK)=PTM1(JL,JK)+PTTENT(JL,JK)*PTSPHY
    ZQP1(JL,JK)=PQM1(JL,JK)+PTTENQ(JL,JK)*PTSPHY
    ZQSAT(JL,JK)=ZQP1(JL,JK)
  ENDDO
ENDDO

IF ( LDSLPHY ) THEN
  IF(KTRAC>0 .AND. LMFTRAC) THEN
    DO JN=1,KTRAC
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
! attention: transport is for positive definite quantities only
!          ZCP1(JL,JK,JN)=MAX(0._JPRB,PCM1(JL,JK,JN))
          ZCP1(JL,JK,JN)=PCM1(JL,JK,JN)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF(KTRAC>0 .AND. LMFTRAC) THEN
    DO JN=1,KTRAC
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
! attention: transport is for positive definit quantities only
!          ZCP1(JL,JK,JN)=MAX(0._JPRB,PCM1(JL,JK,JN)+PTENC(JL,JK,JN)*PTSPHY)
          ZCP1(JL,JK,JN)=PCM1(JL,JK,JN)+PTENC(JL,JK,JN)*PTSPHY
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ENDIF

IFLAG=1
IF (LDRAIN1D) THEN
  CALL SATUR_1D (KIDIA , KFDIA , KLON  , NJKT22 , KLEV,&
   & PAP   , ZTP1  , ZQSAT )  
ELSE
  CALL SATUR (KIDIA , KFDIA , KLON  , NJKT22 , KLEV, YDML_PHY_SLIN%YREPHLI%LPHYLIN, &
   & PAP   , ZTP1  , ZQSAT , IFLAG  ) 
ENDIF 

DO JL=KIDIA,KFDIA
  ZRAIN(JL)=0.0_JPRB
  PCAPE(JL)=0.0_JPRB
  LDCUM(JL)=.FALSE.
ENDDO

!-----------------------------------------------------------------------

!*    2.     CALL 'CUMASTR2'(MASTER-ROUTINE FOR SIMPLIFIED 
!*                           CUMULUS PARAMETERIZATION) 
!*           ---------------------------------------------- 

CALL CUMASTRN2 &
 & (  YDML_PHY_SLIN,   YDML_PHY_EC,  &
 & KIDIA,    KFDIA,    KLON,    KLEV,&
 & LDLAND,   LDRAIN1D, PTSPHY,&
 & ZTP1,     ZQP1,     ZUP1,     ZVP1,&
 & PVERVEL,  ZQSAT,    PQHFL,    PAHFS,&
 & PAP,      PAPH,     PGEO,     PGEOH, PGAW,&
 & PTENT,    PTENQ,    PTENU,    PTENV,&
 & LDCUM,    KTYPE,    KCBOT,    KCTOP,&
 & KBOTSC,   LDSC,&
 & ZTU,      ZQU,      PLU,      PLUDE,&
 & ZENTHD,   PFPLCL,   PFPLCN,   ZRAIN,&
 & PMFU,     PMFD,&
 & PMFUDE_RATE,        PMFDDE_RATE,    PCAPE,&
 & KTRAC,    ZCP1,     PTENC,    PSCAV )  

!----------------------------------------------------------------------
!-----------------------------------------------------------------------

!*    3.0       CALL 'CUCCDIA' TO UPDATE CLOUD PARAMETERS FOR RADIATION
!               -------------------------------------------------------

CALL CUCCDIA &
 & (  YDERAD,YDML_PHY_SLIN%YREPHLI,YDML_PHY_EC%YREPHY, &
 & KIDIA,    KFDIA,    KLON,   KLEV,&
 & NSTEP,    KCBOT,    KCTOP,&
 & LDCUM,    ZQU,      PLU,      PMFU,    ZRAIN,&
 & PARPRC,   KTOPC,    KBASEC                   )  

!----------------------------------------------------------------------
! SECTION 4 IS ONLY REQUIRED IF THE DIAGNOSTIC CLOUD SCHEME IS USED :
! I.E., LEPCLD=.FALSE.

!*    4.0          CALL 'CUSTRAT' FOR PARAMETERIZATION OF PBL-CLOUDS
!                  -------------------------------------------------

LLTEST = (.NOT.LEPCLD.AND..NOT.LENCLD2.AND..NOT.LEPCLD2).OR. &
  & (LPHYLIN.AND..NOT.LENCLD2.AND..NOT.LEPCLD2).AND..NOT.LDRAIN1D

IF (LLTEST) THEN

  CALL CUSTRAT &
   & (  YDML_PHY_SLIN%YREPHLI, &
   & KIDIA,    KFDIA,    KLON,   KLEV,&
   & LDCUM,    PTSPHY,   PVDIFTS,&
   & PAP,      PAPH,     PGEO,&
   & ZTP1,     ZQP1,     ZQSAT,   ZENTHS,&
   & PTENT,    PTENQ)  

ENDIF

!---------------------------------------------------------------------

! Extraction of tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENQ(JL,JK)=PTENQ(JL,JK)-PTTENQ(JL,JK)
    PTENT(JL,JK)=PTENT(JL,JK)-PTTENT(JL,JK)
    PTENV(JL,JK)=PTENV(JL,JK)-PTTENV(JL,JK)
    PTENU(JL,JK)=PTENU(JL,JK)-PTTENU(JL,JK)
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!---------------------------------------------------------------------

!*    5.           FLUX COMPUTATIONS
!                  -----------------

DO JL=KIDIA,KFDIA
  PDIFCQ(JL,1)=0.0_JPRB
  PDIFCS(JL,1)=0.0_JPRB
  PSTRCU(JL,1)=0.0_JPRB
  PSTRCV(JL,1)=0.0_JPRB
  PFCQLF(JL,1)=0.0_JPRB
  PFCQIF(JL,1)=0.0_JPRB
  PFHPCL(JL,1)=0.0_JPRB
  PFHPCN(JL,1)=0.0_JPRB
  PDIFCS(JL,1)=0.0_JPRB
  PDIFCQ(JL,1)=0.0_JPRB
  ZCONDFLL(JL)=0.0_JPRB
  ZCONDFLN(JL)=0.0_JPRB
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

    ZGDPH   = -RG/(      PAPHM1(JL,JK+1)-PAPHM1(JL,JK) )
!...increment in dry static energy is converted to flux of d.s.e.
    ZCP=RCPD*(1+RVTMP2*ZQP1(JL,JK))
    PDIFCS(JL,JK+1)=PTENT(JL,JK)/ZGDPH*ZCP+ PDIFCS(JL,JK)
!...increment in Q
    PDIFCQ(JL,JK+1)=PTENQ(JL,JK)/ZGDPH + PDIFCQ(JL,JK)
!...increments in U,V
    PSTRCU(JL,JK+1)=PTENU(JL,JK)/ZGDPH + PSTRCU(JL,JK)
    PSTRCV(JL,JK+1)=PTENV(JL,JK)/ZGDPH + PSTRCV(JL,JK)

  ENDDO
ENDDO

DO JK=1,KLEV
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA
!... enthalpy flux due to precipitations

    PFHPCL(JL,JK+1)=-PFPLCL(JL,JK+1)*RLVTT
    PFHPCN(JL,JK+1)=-PFPLCN(JL,JK+1)*RLSTT
    IF (LPHYLIN .OR. LDRAIN1D) THEN
      ZALFAW=MIN(1.0_JPRB,0.545_JPRB*(TANH(0.17_JPRB*(ZTP1(JL,JK)-RLPTRC))+1.0_JPRB))
    ELSE
      ZALFAW=FOEALFCU(ZTP1(JL,JK))
    ENDIF
    ZCONDFLL(JL)=ZCONDFLL(JL) + ZALFAW     *PLUDE(JL,JK)
    ZCONDFLN(JL)=ZCONDFLN(JL) + (1.0_JPRB-ZALFAW)*PLUDE(JL,JK)
    PDIFCS(JL,JK+1)=PDIFCS(JL,JK+1)-PFHPCL(JL,JK+1)-PFHPCN(JL,JK+&
     & 1)&
     & + RLVTT * ZCONDFLL(JL) + RLSTT * ZCONDFLN(JL)  
    PDIFCQ(JL,JK+1)=PDIFCQ(JL,JK+1)-PFPLCL(JL,JK+1)-PFPLCN(JL,JK+&
     & 1)&
     & - ZCONDFLL(JL) - ZCONDFLN(JL)  
    PFCQLF(JL,JK+1)=ZCONDFLL(JL)
    PFCQIF(JL,JK+1)=ZCONDFLN(JL)
  ENDDO
ENDDO

!---------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUCALLN2',1,ZHOOK_HANDLE)
END SUBROUTINE CUCALLN2
