! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUCALLN &
 & (PPLDARE, PPLRG,    YDTHF,    YDCST,   YDMODEL, &
 & KIDIA,    KFDIA,    KLON,     KSMAX,   KLEV, &
 & LDLAND, LDSLPHY,&
 & PTSPHY,PVDIFTS,&
 & PTM1,     PQM1,     PUM1,     PVM1,    PLITOT,&
 & PVERVEL,  PQHFL,    PAHFS,   PAPHM1,&
 & PAP,      PAPH,     PGEO,     PGEOH, PGAW,& 
 & PCUCONVCA,PGP2DSPP, &
 & TENDENCY_CML,TENDENCY_LOC, TENDENCY_DYN,&
 & PTENDENCY_VD9,&
 & PARPRC,&
 & KTOPC,    KBASEC,   KTYPE,&
 & KCBOT,    KCTOP,    KBOTSC,   LDCUM,   LDSC,&
 & KCBOT_LIG, KCTOP_LIG, LDCUM_LIG,&
 & LDSHCV,&
 & PLCRIT_AER,&     
 & PLU,      PLUDE,    PLUDELI,  PSNDE,    PMFU,    PMFD,&
 & PDIFCQ,   PDIFCS,   PFHPCL,   PFHPCN,&
 & PFPLCL,   PFPLCN,   PLRAIN,   PRSUD,&
 & PSTRCU,   PSTRCV,   PFCQLF,   PFCQIF,&
 & PMFUDE_RATE ,       PMFDDE_RATE ,      PCAPE  ,PWMEAN, PVDISCU, PDISS,  &
 & KTRAC,    PCM1,     PTENC,    PSCAV )  

!          *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                     *CUMASTR* (CUMULUS PARAMETERIZATION)
!                     *CUCCDIA* (CUMULUS CLOUDS FOR RADIATION)
!                     *CUSTRAT* (PBL_STRATOCUMULUS)

!**   PURPOSE.
!     --------

!          *CUCALL* - INTERFACE FOR *CUMASTR*,*CUCCDIA* AND *CUSTRAT*:
!                     PROVIDES INPUT FOR CUMASTR, CUCCDIA AND CUSTRAT.
!                     RECEIVES UPDATED TENDENCIES, PRECIPITATION
!                     AND CONVECTIVE CLOUD PARAMETERS FOR RADIATION.

!**   INTERFACE.
!     ----------

!          *CUCALL* IS CALLED FROM *CALLPAR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KSMAX*        HORIZONTAL SPECTRAL TRUNCATION
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
!    *PLITOT*       GRID MEAN LIQUID WATER + ICE CONTENT          KG/KG
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
!    *PLCRIT_AER*   CRITICAL LIQUID MMR FOR AUTOCONVERSION PROCESS KG/KG
!    *PCUCONVCA*    COUPLING TO CELL AUTOMATON OR 2D ADVECT FIELD ( )
!    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)

!    UPDATED PARAMETERS (REAL):

!     tendency_cml    cumulative tendency used for final output
!     tendency_loc    local tendency from convection
!     tendency_dyn    tendency from dynamics
!     ptendency_vd9   tendency from vertical diffusion stored from previous timestep

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
!    *KCBOT_LIG*    CLOUD BASE LEVEL (FOR LIGHTNING PARAM)
!    *KCTOP_LIG*    CLOUD TOP LEVEL (FOR LIGHTNING PARAM)

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDCUM_LIG*    FLAG: .TRUE. FOR CONVECTIVE POINTS (FOR LIGHTNING PARAM)
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (REAL):

!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
!    *PLUDE*        DETRAINED TOTAL CONDENSATE                  KG/(M2*S)
!    *PLUDELI*      DETRAINED LIQUID, ICE, VAPOR, T             KG/(M2*S)
!    *PSNDE*        DETRAINED SNOW/RAIN                         KG/(M2*S)
!    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                         KG/(M2*S)
!    *PDIFCQ*       CONVECTIVE FLUX OF SPECIFIC HUMIDITY        KG/(M2*S)
!    *PDIFCS*       CONVECTIVE FLUX OF HEAT                      J/(M2*S)
!    *PFHPCL*       ENTHALPY FLUX DUE TO CONVECTIVE RAIN         J/(M2*S)
!    *PFHPCN*       ENTHALPY FLUX DUE TO CONVECTIVE SNOW         J/(M2*S)
!    *PFPLCL*       CONVECTIVE RAIN FLUX                        KG/(M2*S)
!    *PFPLCN*       CONVECTIVE SNOW FLUX                        KG/(M2*S)
!    *PLRAIN*       RAIN+SNOW CONTENT IN UPDRAFTS                 KG/KG
!    *PRSUD*        CONVECT MEAN RAIN AND SNOW CONTENT IN COLUMN  KG/KG
!    *PSTRTU*       CONVECTIVE FLUX OF U-MOMEMTUM         (M/S)*KG/(M2*S)
!    *PSTRTV*       CONVECTIVE FLUX OF V-MOMEMTUM         (M/S)*KG/(M2*S)
!    *PFCQLF*       CONVECTIVE CONDENSATION FLUX LIQUID         KG/(M2*S)
!    *PFCQIF*       CONVECTIVE CONDENSATION FLUX ICE            KG/(M2*S)

!               Diagnostics:
!    *PMFUDE_RATE* UPDRAFT DETRAINMENT RATE                     KG/(M3*S)
!    *PMFDDE_RATE* DOWNDRAFT DETRAINMENT RATE                   KG/(M3*S)
!    *PCAPE*       MAXIMUM pseudoadiabat CAPE                   (J/KG) 
!    *PWMEAN*      VERTICALLY AVERAGED UPDRAUGHT VELOCITY        M/S
!    *PVDISCU*     VERT INTEGRATED CONVECTIVE DISSIPATION RATE   W/m2
!    *PDISS*       ESTIMATE FOR DISSIPATION BY DEEP CONVCTION    m2/s3^(1/3)

!     EXTERNALS.
!     ----------

!          CUMASTR
!          CUCCDIA
!          CUSTRAT

!     AUTHOR.
!     -------
!      M.TIEDTKE      E.C.M.W.F.     12/1989

!     MODIFICATIONS.
!     --------------
!      15-10-01 : FULLIMP mods D.Salmond 
!      11-02-04 : Enable tracer transport P. Bechtold
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!      F.Vana        18-May-2012 cleaning
!      P.Lopez       08-Dec-2020 Added separate fields for lightning parameterization.
!      R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
!      F.Vana&P.Bechtold (Jan 2023): new seq. physics order
!----------------------------------------------------------------------

USE TYPE_MODEL                   , ONLY : MODEL
USE PARKIND1                     , ONLY : JPIM     ,JPRB
USE YOMHOOK                      , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST                       , ONLY : TCST  
USE YOETHF                       , ONLY : TTHF  
USE YOMCT3                       , ONLY : NSTEP
USE YOMCAPE                      , ONLY : LMCAPEA
USE YOMPHYDER                    , ONLY : STATE_TYPE

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PPLDARE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPLRG
TYPE(TTHF)        ,INTENT(IN)    :: YDTHF
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON) 
LOGICAL           ,INTENT(IN)    :: LDSLPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLCRIT_AER(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLITOT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUCONVCA(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCAV(KTRAC) 
TYPE (STATE_TYPE) ,INTENT(IN)    :: TENDENCY_CML   ! cumulative tendency used for final output
TYPE (STATE_TYPE) ,INTENT(INOUT) :: TENDENCY_LOC   ! local tendency from convection
TYPE (STATE_TYPE) ,INTENT(IN)    :: TENDENCY_DYN   ! dynamics tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_VD9(KLON,KLEV,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND_VD)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PARPRC(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KTOPC(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBASEC(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KTYPE(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCBOT_LIG(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCTOP_LIG(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDCUM(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDCUM_LIG(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON) 
LOGICAL           ,INTENT(IN)    :: LDSHCV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLUDE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLUDELI(KLON,KLEV,4) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSNDE(KLON,KLEV,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCQ(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLRAIN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSUD(KLON,KLEV,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCV(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQLF(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQIF(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFUDE_RATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFDDE_RATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWMEAN(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDISCU(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISS(KLON,KLEV)

REAL(KIND=JPRB) ::     ZTP1(KLON,KLEV),        ZQP1(KLON,KLEV),&
 & ZUP1(KLON,KLEV),        ZVP1(KLON,KLEV),&
 & ZTU(KLON,KLEV),         ZQU(KLON,KLEV),&
 & ZQSAT(KLON,KLEV),       ZRAIN(KLON)         

REAL(KIND=JPRB) :: ZCP1(KLON,KLEV,KTRAC)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZENTHD(KLON,KLEV),ZENTHS(KLON,KLEV),ZDX(KLON)  
REAL(KIND=JPRB) :: ZCONDFLL(KLON),ZCONDFLN(KLON)

INTEGER(KIND=JPIM) :: IFLAG, JK, JL, JN

REAL(KIND=JPRB) :: ZCP, ZGDPH, ZEPS

LOGICAL :: LLTEST, LLTDKMF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuccdia.intfb.h"
#include "cumastrn.intfb.h"
#include "custrat.intfb.h"
#include "satur.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "fcttrm.func.h"

!-----------------------------------------------------------------------

!*    0.1          STORE T,Q,X,Y TENDENCIES FOR FLUX COMPUTATIONS
!*                 ----------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCALLN',0,ZHOOK_HANDLE)
ASSOCIATE(YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDML_PHY_SLIN=>YDMODEL%YRML_PHY_SLIN,&
 & YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI,YDML_PHY_EC=>YDMODEL%YRML_PHY_EC,&
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM,&
 & YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY,YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG)
ASSOCIATE(NJKT2=>YDML_PHY_EC%YRECUMF%NJKT2, &
 & LPHYLIN=>YDML_PHY_SLIN%YREPHLI%LPHYLIN, RLPTRC=>YDML_PHY_SLIN%YREPHLI%RLPTRC, &
 & RCPD=>YDCST%RCPD, RG=>YDCST%RG, RLSTT=>YDCST%RLSTT, RLVTT=>YDCST%RLVTT, &
 & RVTMP2=>YDTHF%RVTMP2, &
 & LEPCLD=>YDML_PHY_EC%YREPHY%LEPCLD, LMFTRAC=>YDML_PHY_EC%YREPHY%LMFTRAC, &
 & LENCLD2=>YDML_PHY_SLIN%YRPHNC%LENCLD2,LEPCLD2=>YDML_PHY_SLIN%YRPHNC%LEPCLD2, &
 & MT_SAVTEND=>YDSLPHY%MT_SAVTEND, MU_SAVTEND=>YDSLPHY%MU_SAVTEND, &
 & MV_SAVTEND=>YDSLPHY%MV_SAVTEND, MQ_SAVTEND=>YDSLPHY%MQ_SAVTEND, &
 & REPQMI=>YDML_PHY_EC%YRECND%REPQMI)

! Setup of tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    TENDENCY_LOC%Q(JL,JK)=TENDENCY_CML%Q(JL,JK)+PTENDENCY_VD9(JL,JK,MQ_SAVTEND)
    TENDENCY_LOC%T(JL,JK)=TENDENCY_CML%T(JL,JK)+PTENDENCY_VD9(JL,JK,MT_SAVTEND)
    TENDENCY_LOC%V(JL,JK)=TENDENCY_CML%V(JL,JK)+PTENDENCY_VD9(JL,JK,MV_SAVTEND)
    TENDENCY_LOC%U(JL,JK)=TENDENCY_CML%U(JL,JK)+PTENDENCY_VD9(JL,JK,MU_SAVTEND)
    ZENTHD(JL,JK)=0.0_JPRB
    ZENTHS(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    1.           UPDATE PROGN. VALUES DEPENDING ON 
!                  TIME INTEGRATION AND CALCULATE QS
!*                 ---------------------------------
!-------------------------------------------------
! Protection in division. ZEPS is choosen as 1000 times lower
! than the actual values of qv (1.e-7 kg/kg) met in the high atmosphere.
! This ZEPS value, used in multiplication only, is relevant for both simple and
! double precision purpose.
!-------------------------------------------------
ZEPS=1.0E-10_JPRB
!
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZUP1(JL,JK)=PUM1(JL,JK)+TENDENCY_CML%U(JL,JK)*PTSPHY
    ZVP1(JL,JK)=PVM1(JL,JK)+TENDENCY_CML%V(JL,JK)*PTSPHY
    ZTP1(JL,JK)=PTM1(JL,JK)+TENDENCY_CML%T(JL,JK)*PTSPHY
    ZQP1(JL,JK)=MAX(REPQMI,PQM1(JL,JK)+TENDENCY_CML%Q(JL,JK)*PTSPHY)
    ZQSAT(JL,JK)=ZQP1(JL,JK)
    IF(ZQSAT(JL,JK)<=0.0_JPRB) ZQSAT(JL,JK)=ZEPS
  ENDDO
ENDDO

! NOTE: Following has not much to do with LSLPHY. It works for NL code
!       but not really in TL/AD. Perhaps there should be a better key
!       to uniformely control following.
IF ( LDSLPHY ) THEN
  IF(KTRAC>0 .AND. LMFTRAC) THEN
    DO JN=1,KTRAC
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
   ! attention: transport is for positive definit quantities only
        ! ZCP1(JL,JK,JN)=MAX(0._JPRB,PCM1(JL,JK,JN))
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
        ! ZCP1(JL,JK,JN)=MAX(0._JPRB,PCM1(JL,JK,JN)+PTENC(JL,JK,JN)*PTSPHY)
          ZCP1(JL,JK,JN)=PCM1(JL,JK,JN)+PTENC(JL,JK,JN)*PTSPHY
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ENDIF

IFLAG=1
CALL SATUR (YDTHF, YDCST, KIDIA , KFDIA , KLON  , NJKT2 , KLEV,&
 & LPHYLIN, &
 & PAP   , ZTP1  , ZQSAT , IFLAG  )  

!-----------------------------------------------------------------------

!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION) 
!*           ----------------------------------------------------------- 


LLTDKMF = .FALSE.
CALL CUMASTRN &
 & (PPLDARE, PPLRG,    YDTHF,   YDCST,    YDML_PHY_SLIN,   YDML_PHY_EC,   YGFL, YDCHEM, YDSPP_CONFIG, &
 & KIDIA,    KFDIA,    KLON,    KLEV, ZDX, LLTDKMF, &
 & LMCAPEA,  LDLAND,   PTSPHY,&
 & ZTP1,     ZQP1,     ZUP1,     ZVP1,     PLITOT,&
 & PVERVEL,  ZQSAT,    PQHFL,    PAHFS,&
 & PAP,      PAPH,     PGEO,     PGEOH,  PGAW,&
 & PCUCONVCA,PGP2DSPP, &
 & TENDENCY_LOC%T,    TENDENCY_LOC%Q,    TENDENCY_LOC%U,    TENDENCY_LOC%V,&
 & TENDENCY_DYN%T,    TENDENCY_DYN%Q, &
 & LDCUM,    KTYPE,    KCBOT,    KCTOP,&
 & LDCUM_LIG, KCBOT_LIG, KCTOP_LIG,&
 & KBOTSC,   LDSC,&
 & LDSHCV,&
 & PLCRIT_AER,&
 & ZTU,      ZQU,      PLU,      PLUDE,  PLUDELI,   PSNDE,&
 & ZENTHD,   PFPLCL,   PFPLCN,   ZRAIN,  PLRAIN,    PRSUD,&
 & PMFU,     PMFD,&
 & PMFUDE_RATE,        PMFDDE_RATE,    PCAPE,  PWMEAN, PVDISCU, PDISS, &
 & KTRAC,    ZCP1,     PTENC,    PSCAV)  

!----------------------------------------------------------------------

!*    3.0       CALL 'CUCCDIA' TO UPDATE CLOUD PARAMETERS FOR RADIATION
!               -------------------------------------------------------

CALL CUCCDIA &
 & (YDERAD,  YDEPHLI,  YDEPHY,&
 & KIDIA,    KFDIA,    KLON,   KLEV,&
 & NSTEP,    KCBOT,    KCTOP,&
 & LDCUM,    ZQU,      PLU,      PMFU,    ZRAIN,&
 & PARPRC,   KTOPC,    KBASEC                   )  

!----------------------------------------------------------------------
! SECTION 4 IS ONLY REQUIRED IF THE DIAGNOSTIC CLOUD SCHEME IS USED :
! I.E., LEPCLD=.FALSE.

LLTEST = (.NOT.LEPCLD.AND..NOT.LENCLD2.AND..NOT.LEPCLD2).OR. &
 & (LPHYLIN.AND..NOT.LENCLD2.AND..NOT.LEPCLD2)

IF (LLTEST) THEN
!---------------------------------------------------------------------

!*    4.0          CALL 'CUSTRAT' FOR PARAMETERIZATION OF PBL-CLOUDS
!                  -------------------------------------------------

  CALL CUSTRAT &
   & (  YDEPHLI,&
   & KIDIA,    KFDIA,    KLON,   KLEV,&
   & LDCUM,    PTSPHY,   PVDIFTS,&
   & PAP,      PAPH,     PGEO,&
   & ZTP1,     ZQP1,     ZQSAT,   ZENTHS,&
   & TENDENCY_LOC%T,    TENDENCY_LOC%Q)  

ENDIF

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

    TENDENCY_LOC%Q(JL,JK)=TENDENCY_LOC%Q(JL,JK)-PTENDENCY_VD9(JL,JK,MQ_SAVTEND)
    TENDENCY_LOC%T(JL,JK)=TENDENCY_LOC%T(JL,JK)-PTENDENCY_VD9(JL,JK,MT_SAVTEND)
    TENDENCY_LOC%V(JL,JK)=TENDENCY_LOC%V(JL,JK)-PTENDENCY_VD9(JL,JK,MV_SAVTEND)
    TENDENCY_LOC%U(JL,JK)=TENDENCY_LOC%U(JL,JK)-PTENDENCY_VD9(JL,JK,MU_SAVTEND)

    ZGDPH   = -RG/(      PAPHM1(JL,JK+1)-PAPHM1(JL,JK) )
!...increment in dry static energy is converted to flux of d.s.e.
    ZCP=RCPD*(1+RVTMP2*ZQP1(JL,JK))
    PDIFCS(JL,JK+1)=(TENDENCY_LOC%T(JL,JK)-TENDENCY_CML%T(JL,JK))/ZGDPH*ZCP+ PDIFCS(JL,JK)
!...increments in U,V
    PSTRCU(JL,JK+1)=(TENDENCY_LOC%U(JL,JK)-TENDENCY_CML%U(JL,JK))/ZGDPH + PSTRCU(JL,JK)
    PSTRCV(JL,JK+1)=(TENDENCY_LOC%V(JL,JK)-TENDENCY_CML%V(JL,JK))/ZGDPH + PSTRCV(JL,JK)

!...increment in Q
    PDIFCQ(JL,JK+1)=(TENDENCY_LOC%Q(JL,JK)-TENDENCY_CML%Q(JL,JK))/ZGDPH + PDIFCQ(JL,JK)

  ENDDO
ENDDO

DO JK=1,KLEV
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA
!... enthalpy flux due to precipitations

    PFHPCL(JL,JK+1)=-PFPLCL(JL,JK+1)*RLVTT
    PFHPCN(JL,JK+1)=-PFPLCN(JL,JK+1)*RLSTT
    ZCONDFLL(JL)=ZCONDFLL(JL) + PLUDELI(JL,JK,1)
    ZCONDFLN(JL)=ZCONDFLN(JL) + PLUDELI(JL,JK,2)
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

! Extraction of tendencies relevant to the convection
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    TENDENCY_LOC%Q(JL,JK)=TENDENCY_LOC%Q(JL,JK)-TENDENCY_CML%Q(JL,JK)
    TENDENCY_LOC%T(JL,JK)=TENDENCY_LOC%T(JL,JK)-TENDENCY_CML%T(JL,JK)
    TENDENCY_LOC%V(JL,JK)=TENDENCY_LOC%V(JL,JK)-TENDENCY_CML%V(JL,JK)
    TENDENCY_LOC%U(JL,JK)=TENDENCY_LOC%U(JL,JK)-TENDENCY_CML%U(JL,JK)
  ENDDO
ENDDO
!---------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUCALLN',1,ZHOOK_HANDLE)
END SUBROUTINE CUCALLN
