! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_DIAG1 &
 &( YDERAD,YDECLDP,YDML_GCONF,YDPHY2,KIDIA, KFDIA, KLON , KTDIA, KLEV , KLEVX, KFLDX, KIEXT3D, &
 &  PAERO, PAPH , PAP  , PA   , PI   , PL   , & 
 &  PLSM , PQ   , PQSAT, PT   , PWND , &
 &  PEXTRA &
 &)

!**** *AER_DIAG1* - ADDITIONAL AER/CLOUD DIAGNOSTICS FOR MACC

!**   INTERFACE.
!     ----------

!     EXPLICIT ARGUMENTS :
!    ---------------------

! KLEVX   : NUMBER OF LEVELS IN EXTRA MULTI-LEVEL DIAGNOSTIC FIELDS
! KFLDX   : NUMBER OF VARIABLES IN EXTRA MULTI-LEVEL DIAGNOSTIC FIELDS

! PAERO : (KLON,KLEV,NAERO) ; PROGNOSTIC AEROSOL MIXING RATIO
! PAPH  : (KLON,KLEV+1)     ; HALF-LEVEL PRESSURE
! PAP   : (KLON,KLEV)       ; FULL LEVEL PRESSURE
! PA    : (KLON,KLEV)       ; CLOUD FRACTION
! PI    : (KLON,KLEV)       ; CLOUD ICE MIXING RATIO
! PL    : (KLON,KLEV)       ; CLOUD LIQUID WATER MIXING RATIO
! PLSM  : (KLON)            ; LAND-SEA MASK
! PQ    : (KLON,KLEV)       ; SPECIFIC HUMIDITY
! PQSAT : (KLON,KLEV)       ; SATURATION HUMIDITY
! PT    : (KLON,KLEV)       ; TEMPERATURE
! PWND  : (KLON)            ; 10-M WIND

! PEXTRA:                   ; EXTRA DIAGNOSTIC FIELD

   
!     AUTHORS.
!     --------
!        J.-J. MORCRETTE         *ECMWF*
!                       USING MATERIAL FROM A.TOMPKINS' *AER_CLCLD*

!     Modifications:
!     --------------
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMPHY2  , ONLY : TPHY2
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD, RG, RTT, RPI
USE YOMCT3   , ONLY : NSTEP
USE YOECLDP  , ONLY : TECLDP
USE YOERAD   , ONLY : TERAD
USE YOERDU   , ONLY : REPLOG, REPSCA
USE YOETHF   , ONLY : RTICE
USE YOMLUN   , ONLY : NULOUT

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TPHY2)       ,INTENT(INOUT) :: YDPHY2
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX
INTEGER(KIND=JPIM),INTENT(IN)    :: KIEXT3D

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV,YDML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWND(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX) 

!     -----------------------------------------------------------------

INTEGER(KIND=JPIM) :: IRADIP, IRADLP
INTEGER(KIND=JPIM) :: JAERSS, JAERDD , JAEROM, JAERBC, JAERSU
INTEGER(KIND=JPIM) :: JAER, JL, JK

LOGICAL :: LLPRINT, LLCLOUD(KLEV)

REAL(KIND=JPRB) :: ZCLFR(KLON,KLEV), ZFIWP(KLON,KLEV), ZFLWP(KLON,KLEV) 
REAL(KIND=JPRB) :: ZPS(KLON)       , ZQIWP(KLON,KLEV), ZQLWP(KLON,KLEV)
REAL(KIND=JPRB) :: ZCCNL(KLON)     , ZCCNO(KLON)     , ZCCN0(KLON)
REAL(KIND=JPRB) :: ZQ(KLON)        , ZEXPL(KLON)     , ZEXPO(KLON)
REAL(KIND=JPRB) :: ZRHO(KLON,KLEV) , ZRE_LIQ(KLON,KLEV)
REAL(KIND=JPRB) :: ZDESR(KLON)     , ZRADIP(KLON)    , ZRADLP(KLON)
REAL(KIND=JPRB) :: ZMAERMN(5)      , ZMAER(KLON,KLEV,5)
REAL(KIND=JPRB) :: ZAERO(KLON,KLEV,YDML_GCONF%YGFL%NACTAERO)         , ZCCN(KLON,KLEV)
REAL(KIND=JPRB) :: ZICENUCLEI(KLON,KLEV), ZNICE(KLON,KLEV), ZRE_ICE(KLON,KLEV)

! ZIWC  : (KLON,KLEV)       ; CLOUD ICE CONTENT
! ZWLC  : (KLON,KLEV)       ; CLOUD LIQUID WATER CONTENT
! ZDDE  : (KLON,KLEV)       ; ICE PARTICLE DIMENSION (DIAGNOSTIC)
! ZPDE  : (KLON,KLEV)       ; ICE PARTICLE DIMENSION (PROGNOSTIC)
! ZDRE  : (KLON,KLEV)       ; LIQUID WATER CLOUD EFFECTIVE RADIUS (DIAGNOSTIC)
! ZPRE  : (KLON,KLEV)       ; LIQUID WATER CLOUD EFFECTIVE RADIUS (PROGNOSTIC)
! ZCCND : (KLON,KLEV)       ; CLOUD CONDENSATION NUCLEI (DIAGNOSTIC)
! ZCCNP : (KLON,KLEV)       ; CLOUD CONDENSATION NUCLEI (PROGNOSTIC)
! ZDNC1 : (KLON,KLEV)       ; CLOUD DROPLET NUMBER CONCENTRATION (DIAGNOSTIC)
! ZDNC2 : (KLON,KLEV)       ; CLOUD DROPLET NUMBER CONCENTRATION (PROGNOSTIC)
REAL(KIND=JPRB) :: ZIWC(KLON,KLEV)
REAL(KIND=JPRB) :: ZLWC(KLON,KLEV)
REAL(KIND=JPRB) :: ZDDE(KLON,KLEV) , ZPDE(KLON,KLEV)
REAL(KIND=JPRB) :: ZDRE(KLON,KLEV) , ZPRE(KLON,KLEV)
REAL(KIND=JPRB) :: ZCCND(KLON,KLEV), ZCCNP(KLON,KLEV)
REAL(KIND=JPRB) :: ZDNC1(KLON,KLEV), ZDNC2(KLON,KLEV)

REAL(KIND=JPRB) :: ZDP  , ZDPOG , ZIWGKG, ZLWGKG, ZPODT, ZREFDE, ZDEFRE
REAL(KIND=JPRB) :: ZASEA, ZALND , ZD    , ZNTOT , ZNUM , ZDEN  , ZTEMP
REAL(KIND=JPRB) :: ZRG  , ZTEMPC, ZTCELS, ZFSR  , ZAIWC, ZBIWC , ZEXPN
REAL(KIND=JPRB) :: ZRELVNT, ZVOLD
REAL(KIND=JPRB) :: ZEPSEC, ZS0, ZSVP, ZSCRITHOMO, ZNICEHOMO
REAL(KIND=JPRB) :: ZCLD, ZCLDMAX, ZNCRIT_GIERENS, ZNCRIT_REN, ZRHO_ICE, ZRHO_LIQ, ZWTOT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------------

#include "update_fields.intfb.h"

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AER_DIAG1',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YDML_GCONF%YGFL%NACTAERO, &
 & RCCNOM=>YDECLDP%RCCNOM, RCCNSS=>YDECLDP%RCCNSS, RCCNSU=>YDECLDP%RCCNSU, &
 & RNICE=>YDECLDP%RNICE, &
 & LCCNL=>YDERAD%LCCNL, LCCNO=>YDERAD%LCCNO, NRADIP=>YDERAD%NRADIP, &
 & NRADLP=>YDERAD%NRADLP, RCCNLND=>YDERAD%RCCNLND, RCCNSEA=>YDERAD%RCCNSEA, &
 & RRE2DE=>YDERAD%RRE2DE, &
 & NSTART=>YDML_GCONF%YRRIP%NSTART)
!     -----------------------------------------------------------------

!*         0.     PREPARATORY CALCULATIONS
!                 ------------------------

! ZWTOT is not set !!!
ZWTOT=0.0_JPRB
LLPRINT=.FALSE.
IF (NSTEP < NSTART + 1) THEN
  LLPRINT=.TRUE.
ENDIF
LLPRINT=.FALSE.
IF (LLPRINT) WRITE(NULOUT,FMT='(" NRADLP NRADIP ",2I4,E12.5,1X,2L3)') NRADLP,NRADIP,RRE2DE,LCCNL,LCCNO

ZIWC(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZLWC(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZDDE(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZPDE(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZDRE(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZPRE(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZCCND(KIDIA:KFDIA,1:KLEV)=0._JPRB
ZCCNP(KIDIA:KFDIA,1:KLEV)=0._JPRB
ZDNC1(KIDIA:KFDIA,1:KLEV)=0._JPRB
ZDNC2(KIDIA:KFDIA,1:KLEV)=0._JPRB

!-------------------------------------------------------------------------------

!*         0.5    SECURITY CALCULATIONS
!                 ---------------------

ZEPSEC=1.E-10_JPRB
ZRHO_ICE=900._JPRB
ZRHO_LIQ=1000._JPRB
ZCLDMAX=5.E-03_JPRB

ZRG=1.0_JPRB/RG
DO JAER=1,NACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAERO(JL,JK,JAER)=MAX(0._JPRB,PAERO(JL,JK,JAER))
    ENDDO
  ENDDO
ENDDO

!-------------------------------------------------------------------------------

!*         1.     IN-CLOUD ICE AND LIQUID WATER MIXING RATIOS
!                 -------------------------------------------

!-- in-cloud ice and water mixing ratios; cloud ice and liquid contents and paths

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZCLFR(JL,JK)=PA(JL,JK)
    IF (ZCLFR(JL,JK) >= 0.001_JPRB) THEN
      ZTEMP=1.0_JPRB/PA(JL,JK)
      ZQIWP(JL,JK)=MAX(0._JPRB, PI(JL,JK)*ZTEMP)
      ZQLWP(JL,JK)=MAX(0._JPRB, PL(JL,JK)*ZTEMP)
      ZIWGKG=ZQIWP(JL,JK)*1000._JPRB
      ZLWGKG=ZQLWP(JL,JK)*1000._JPRB 
    ELSE
      ZCLFR(JL,JK)=0._JPRB
      ZQIWP(JL,JK)=0._JPRB
      ZQLWP(JL,JK)=0._JPRB
      ZIWGKG=0._JPRB
      ZLWGKG=0._JPRB
    ENDIF
    ZDP=PAPH(JL,JK+1)-PAPH(JL,JK)
    ZDPOG=ZDP * ZRG
    ZFIWP(JL,JK) = ZIWGKG*ZDPOG
    ZFLWP(JL,JK) = ZLWGKG*ZDPOG
    ZPODT = PAP(JL,JK)/(RD*PT(JL,JK))
    ZRHO(JL,JK) = ZPODT
!-- ice and liquid water content in g m-3
    ZIWC(JL,JK) = ZIWGKG*ZPODT
    ZLWC(JL,JK) = ZLWGKG*ZPODT
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  ZPS(JL)=PAPH(JL,KLEV+1)
ENDDO

!-------------------------------------------------------------------------------

!*         2.     DIAGNOSTIC CCN FROM 10-M WIND
!                 -----------------------------

!-- original bounds from Martin et al. (1994) are 375/1500 and 36/280
!-- bounds from the relationship between 10m wind and CCNs
DO JL=KIDIA,KFDIA
!-- over ocean
  IF (PWND(JL) > 30._JPRB) THEN
    ZQ(JL)=327._JPRB  
  ELSEIF (PWND(JL) > 15._JPRB) THEN  
    ZQ(JL)=EXP(0.13_JPRB*PWND(JL)+1.89_JPRB)
  ELSE
    ZQ(JL)=EXP(0.16_JPRB*PWND(JL)+1.44_JPRB)
  ENDIF  
  ZEXPO(JL)=1.2_JPRB+0.5_JPRB*LOG10(ZQ(JL))
  ZCCNO(JL)=10._JPRB**ZEXPO(JL)
  ZCCNO(JL) = MIN(  287._JPRB,MAX( 32._JPRB, ZCCNO(JL) ))

!-- over land
  IF (PWND(JL) <= 15._JPRB) THEN
    ZQ(JL)=EXP(0.16_JPRB*PWND(JL)+1.45_JPRB)
  ELSE
    ZQ(JL)=EXP(0.13_JPRB*PWND(JL)+1.89_JPRB)
  ENDIF
  ZEXPL(JL)=2.21_JPRB+0.3_JPRB*LOG10(ZQ(JL))
  ZCCNL(JL)=10._JPRB**ZEXPL(JL)
  ZCCNL(JL) = MIN( 1743._JPRB,MAX(292._JPRB, ZCCNL(JL) ))
ENDDO

LLCLOUD(1:KLEV)=.FALSE.
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (PA(JL,JK) >= 0.001_JPRB .AND. PL(JL,JK) > 0._JPRB) THEN
      LLCLOUD(JK)=.TRUE.
      ZCCND(JL,JK)=PLSM(JL)*ZCCNL(JL)+(1._JPRB-PLSM(JL))*ZCCNO(JL)
    ENDIF 
  ENDDO
ENDDO

!-------------------------------------------------------------------------------

!*         3.     DIAGNOSTIC EFFECTIVE RADII AND PARTICLE SIZE
!                 --------------------------------------------

ZREFDE = RRE2DE 
ZDEFRE = 1.0_JPRB / ZREFDE
IRADIP=NRADIP
IRADLP=NRADLP


!*         3.1    DIAGNOSTIC EFFECTIVE RADII OF LIQUID WATER CLOUDS
!                 -------------------------------------------------

DO JK=1,KLEV
  DO JL = KIDIA,KFDIA
! --- EFFECTIVE RADIUS FOR WATER, ICE AND RAIN PARTICLES

! very old parametrization as f(pressure)

    IF (IRADLP == 0) THEN
!-- very old parametrization as f(pressure) ERA-15
      ZRADLP(JL)=10.0_JPRB + (100000.0_JPRB-PAP(JL,JK))*3.5_JPRB

    ELSEIF (IRADLP == 1) THEN
! simple distinction between land (10) and ocean (13) Zhang and Rossow
      IF (PLSM(JL) < 0.5_JPRB) THEN
        ZRADLP(JL)=13.0_JPRB
      ELSE
        ZRADLP(JL)=10.0_JPRB
      ENDIF
      
    ELSEIF (IRADLP == 2) THEN
!--  based on Martin et al., 1994, JAS
      IF (PLSM(JL) < 0.5_JPRB) THEN
        IF (LCCNO) THEN
          ZASEA=ZCCNO(JL)
        ELSE  
          ZASEA=RCCNSEA
        ENDIF  
        ZD=0.33_JPRB
        ZNTOT=-1.15E-03_JPRB*ZASEA*ZASEA+0.963_JPRB*ZASEA+5.30_JPRB
      ELSE
        IF (LCCNL) THEN 
          ZALND=ZCCNL(JL)
        ELSE  
          ZALND=RCCNLND
        ENDIF  
        ZD=0.43_JPRB
        ZNTOT=-2.10E-04_JPRB*ZALND*ZALND+0.568_JPRB*ZALND-27.9_JPRB
      ENDIF
      ZNUM=3.0_JPRB*ZLWC(JL,JK)*(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
      ZDEN=4.0_JPRB*RPI*ZNTOT*(1.0_JPRB+ZD*ZD)**3
      ZTEMP=1.0_JPRB/ZDEN
      IF((ZNUM*ZTEMP) > REPLOG)THEN
        ZRADLP(JL)=100.*EXP(0.333*LOG(ZNUM*ZTEMP))
        ZRADLP(JL)=MAX(ZRADLP(JL), 2.0_JPRB)
        ZRADLP(JL)=MIN(ZRADLP(JL),24.0_JPRB)
      ELSE
        ZRADLP(JL)=4.0_JPRB
      ENDIF

    ELSEIF (IRADLP == 3) THEN
!- effective radius of droplets linked to prognostic aerosol 
!  using Menon et al's CCN (see above)
      ZRADLP(JL) = ZRE_LIQ(JL,JK)
      ZRADLP(JL)=MAX(ZRADLP(JL), 2.0_JPRB)
      ZRADLP(JL)=MIN(ZRADLP(JL),24.0_JPRB)
    ENDIF  

!-- effective radius of liquid water clouds (in um)
    IF (PA(JL,JK) > 0.001_JPRB .AND. PL(JL,JK) > 0._JPRB) THEN
      ZDRE(JL,JK)=ZRADLP(JL)
!-- assuming above radius, diagnose the associated cloud droplet number concentration
!-- volume of droplet (in um^3)
    ZVOLD=4._JPRB*RPI/3._JPRB*ZRADLP(JL)**3._JPRB
!-- factor accounts for ZF1 = 1.E-18 (um^3 to m^3)
!                       ZF2 = 1.E-03 (g m-3 to kg m-3)
!                       ZF3 = 1.E+03  kg m-3 (mass of 1 m3 of water)
!                       ZF4 = 1.E+06  (cm-3 to m-3 )
!                       ZFACT = 1./(ZF1*ZF3/(ZF2*ZF4))
      ZCCNP(JL,JK)=ZLWC(JL,JK)/ZVOLD*1.E+06_JPRB
    ENDIF
  ENDDO
  IF (LLPRINT .AND. LLCLOUD(JK)) THEN
    WRITE(NULOUT,9001) JK,(JL,ZCCND(JL,JK),ZCCNP(JL,JK),JL=KIDIA,KFDIA,15)
9001  FORMAT(1X,'aer_diag1 CCND CCNP:',I3,5(2X,I2,2E10.3))
  ENDIF


!*         3.2    DIAGNOSTIC EFFECTIVE PARTICLE DIMENSION OF ICE CLOUDS
!                 -----------------------------------------------------

  DO JL = KIDIA,KFDIA

! diagnosing the ice particle effective radius/diameter

!- ice particle effective radius =f(T) from Liou and Ou (1994)
 
    IF (PT(JL,JK) < RTICE) THEN
      ZTEMPC=PT(JL,JK)-RTT
    ELSE
      ZTEMPC=RTICE-RTT
    ENDIF
    ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
      & 0.0012_JPRB))    
    
    IF (NRADIP == 0) THEN
!-- fixed 40 micron effective radius
      ZRADIP(JL)= 40.0_JPRB
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
      
    ELSEIF (NRADIP == 1) THEN 
!-- old formulation based on Liou & Ou (1994) temperature (40-130microns)    
      ZRADIP(JL)=MAX(ZRADIP(JL), 40.0_JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),130.0_JPRB)
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
      
    ELSEIF (NRADIP == 2) THEN  
!-- formulation following Jakob, Klein modifications to ice content    
      ZRADIP(JL)=MAX(ZRADIP(JL),30.0_JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),60.0_JPRB)
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
 
    ELSEIF (NRADIP == 3  ) THEN
!- ice particle effective radius =f(T,IWC) from Sun and Rikus (1999)
! revised by Sun (2001)
      IF (ZIWC(JL,JK) > 0.0_JPRB ) THEN
        ZTEMPC = PT(JL,JK)-83.15_JPRB
        ZTCELS = PT(JL,JK)-RTT
        ZFSR = 1.2351_JPRB +0.0105_JPRB * ZTCELS
! Sun, 2001 (corrected from Sun & Rikus, 1999)
        ZAIWC = 45.8966_JPRB * ZIWC(JL,JK)**0.2214_JPRB
        ZBIWC = 0.7957_JPRB * ZIWC(JL,JK)**0.2535_JPRB
        ZDESR(JL) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)
!        ZDESR(JL) = MIN ( MAX( ZDESR(JL), RMINICE), 155.0_JPRB)
        ZDESR(JL) = MIN ( MAX( ZDESR(JL),  1._JPRB), 155.0_JPRB)
        ZRADIP(JL)= ZREFDE * ZDESR(JL)
      ELSE
        ZDESR(JL) = 80.0_JPRB
        ZRADIP(JL)= ZREFDE * ZDESR(JL)
      ENDIF  
    ENDIF
    IF (PA(JL,JK) > 0.001_JPRB .AND. PI(JL,JK) > 0._JPRB) THEN
      ZDDE(JL,JK)=ZDESR(JL)
    ENDIF
  ENDDO
  IF (LLPRINT .AND. LLCLOUD(JK)) THEN
    WRITE(NULOUT,9002) JK,(JL,ZDRE(JL,JK),ZDDE(JL,JK),JL=KIDIA,KFDIA,15)
9002  FORMAT(1X,'aer_diag1 ZDRE ZDDE :',I3,5(2X,I2,2E10.3))
  ENDIF
ENDDO

!-------------------------------------------------------------------------------

!*         4.     CCN AND LIQUID WATER CLOUD Re FROM PROGNOSTIC AEROSOLS
!                 ------------------------------------------------------

!-- first aerosol indirect effect: the relevant prognostic aerosols 
!   (OM, SS, SU) are used as cloud condensation nuclei further used
!   to determine the effective radius of droplets in LIQUID water 
!   clouds 

ZCCN(KIDIA:KFDIA,1:KLEV) =0._JPRB
ZRE_LIQ(KIDIA:KFDIA,1:KLEV)=0._JPRB
!-- IF (LAERCCN ) THEN
  IRADLP=3 
  JAERSS=1
  JAEROM=2
  JAERBC=3
  JAERSU=4
  JAERDD=5
  ZMAERMN(JAERSS)=2.12E-10_JPRB*1.E9_JPRB 
  ZMAERMN(JAERDD)=1.01E-09_JPRB*1.E9_JPRB
  ZMAERMN(JAEROM)=3.05E-11_JPRB*1.E9_JPRB
  ZMAERMN(JAERBC)=3.05E-11_JPRB*1.E9_JPRB
  ZMAERMN(JAERSU)=1.02E-09_JPRB*1.E9_JPRB
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
!-- 1st bin of SS
      ZMAER(JL,JK,1) = ZAERO(JL,JK, 1)*ZRHO(JL,JK)*1.E9_JPRB
!-- hydrophilic OM
      ZMAER(JL,JK,2) = ZAERO(JL,JK, 7)*ZRHO(JL,JK)*1.E9_JPRB
!-- hydrophilic BC
      ZMAER(JL,JK,3) = ZAERO(JL,JK, 9)*ZRHO(JL,JK)*1.E9_JPRB
!-- SO4
      ZMAER(JL,JK,4) = ZAERO(JL,JK,11)*ZRHO(JL,JK)*1.E9_JPRB
!-- 1 st bin of DD
      ZMAER(JL,JK,5) = ZAERO(JL,JK, 4)*ZRHO(JL,JK)*1.E9_JPRB

      ZRELVNT=ZMAER(JL,JK,1)+ZMAER(JL,JK,2)+ZMAER(JL,JK,4)

!---------------------------------------------------------------------
! Turn aerosol mass into a CCN Number concentration for warm rain
! From Menon et al, 2002: JAS, 59, 692-713  Eqns 1a, 1b
!---------------------------------------------------------------------

! CCNxx factors: for org.matter 0.13, sea salt 0.05, and sulphate 0.50
      ZEXPN= RCCNOM*LOG10(MAX(ZMAER(JL,JK,JAEROM),REPSCA)) +&
           & RCCNSS*LOG10(MAX(ZMAER(JL,JK,JAERSS),REPSCA)) +&
           & RCCNSU*LOG10(MAX(ZMAER(JL,JK,JAERSU),REPSCA))

!     ZCCN = N in cm**-3
      ZCCN0(JL) = 10.0_JPRB**(2.41 + ZEXPN)

!-- NB: the prognosed CCN is bounded as the diagnostic one
      IF (ZRELVNT > 0._JPRB .AND. PA(JL,JK) > 0.001_JPRB .AND. PL(JL,JK) > 0._JPRB) THEN
        ZCCN(JL,JK) = MIN( 1743._JPRB, MAX( 32._JPRB, ZCCN0(JL) ))
        ZDNC1(JL,JK)= ZCCN(JL,JK)

! number is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]
!-- ZRE_LIQ for subsequent radiative computations (in um once multiplied by 1.E+06)
        ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387E-10_JPRB*ZRHO(JL,JK)*ZQLWP(JL,JK)/ZCCN(JL,JK))**0.333_JPRB
        ZPRE(JL,JK) = ZRE_LIQ(JL,JK)
      ENDIF
    ENDDO
    IF (LLPRINT .AND. LLCLOUD(JK)) THEN
      WRITE(NULOUT,9003) JK,(JL,ZCCN0(JL),ZQLWP(JL,JK),ZRHO(JL,JK),ZDNC1(JL,JK),ZPRE(JL,JK),JL=KIDIA,KFDIA,15)
9003  FORMAT(1X,'aer_diag1 CCN0 LWP RH0 DNC1 ZPRE:',I3,5(2X,I2,5E10.3))
    ENDIF
  ENDDO
!-- ENDIF

!-------------------------------------------------------------------------------

!*         5.     ICE NUCLEI AND ICE WATER CLOUD Re FROM PROGNOSTIC AEROSOLS
!                 ----------------------------------------------------------

LLCLOUD(1:KLEV)=.FALSE.
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (PA(JL,JK) >= 0.001_JPRB .AND. PI(JL,JK) > 0._JPRB) THEN
      LLCLOUD(JK)=.TRUE.
    ENDIF 
  ENDDO
ENDDO


!-- IF (LAERICN ) THEN
!---------------------------------------------------------------------
! Turn aerosol mass into a Ice Number concentration for ice processes
!---------------------------------------------------------------------
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA

!                       0.01_JPRB is "default" value from
! Demott et al., 1994: JAS, 51, 77-90   -->  ISS = 25 %
! Demott et al., 1997: JGR, 102, 19575-19584
! Meyers et al., 1992: JAM, 31, 708-721  --> ISS = 25 %

! Demott et al. Ice SS=55% or Meyers et al. 1992 JAM, ISS=25%
! In a prognostic scheme this will be function of clear sky humidity

! By relating IN to Aerosol mass we are assuming that the mode of the 
! aerosol size distribution lies in the accumulation or coarse mode

! The relationship will implicitly introduce the exponential height
! dependence that Sassen (1992) and K and Curry (1998) explicitly 
! introduced to their parametrizations.
! Sassen, 1992: 
! Khvorostyanov and Curry, 2000: GRL 27, 4081-4084.

      ZICENUCLEI(JL,JK)=0.01_JPRB*&
    &  (ZMAER(JL,JK,JAERSU)+ZMAER(JL,JK,JAERBC)+ZMAER(JL,JK,JAERDD)) &
    & /(ZMAERMN(JAERSU)    +ZMAERMN(JAERBC)    +ZMAERMN(JAERDD))
      ZICENUCLEI(JL,JK)=MAX(ZICENUCLEI(JL,JK),0.0_JPRB)

! T in oC
      ZTEMPC=PT(JL,JK)-RTT

! Re form for ice crystals from Liou and Oort 1994
! used to derive Re(ice) as in Lohmann, 2002: JAS, 59, 647-656. Eqn 5 
      ZNICEHOMO=0.0_JPRB
      ZRE_ICE(JL,JK)=0.5_JPRB*(326.3_JPRB+ZTEMPC* &
        & (12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*0.0012_JPRB)))
      ZRE_ICE(JL,JK)=MAX(ZRE_ICE(JL,JK),0.0_JPRB)

! effect Re to volume mean from S Moss or Lohmann and Kaercher papers 200?
      ZRE_ICE(JL,JK)=(MAX(SQRT(5.113E6_JPRB+2.809E3_JPRB*ZRE_ICE(JL,JK)**3.0_JPRB)-2.261E3_JPRB,0.0_JPRB))**0.333_JPRB
      ZRE_ICE(JL,JK)=MAX(ZRE_ICE(JL,JK),1.0_JPRB)  ! diameter minimum 1.0 microns

! more default values if not applying
      ZNICE(JL,JK)=RNICE ! place as default

      IF (PT(JL,JK)<238._JPRB .AND. PI(JL,JK)>ZEPSEC) THEN
        ZS0=1.3_JPRB
        ZSCRITHOMO=2.349_JPRB-PT(JL,JK)/259.0_JPRB !ren form of Koop 2000 
        ZSVP=MAX(ZEPSEC,PQSAT(JL,JK)*PAP(JL,JK)/0.622_JPRB)

! Klaus Gierens critical ice nuclei: Gierens, 2003: ACP, 3, 437-446.
        ZNCRIT_GIERENS=2.81E11_JPRB*(10.0_JPRB**(4.0_JPRB-0.02_JPRB*PT(JL,JK)))**0.75_JPRB&
       &*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB/&
       &(PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB*(ZSCRITHOMO-ZS0)**0.75_JPRB)
        ZNCRIT_GIERENS=ZNCRIT_GIERENS/1.E6_JPRB ! cm**-3

! Ren and Mackensie, 2005: QJRMS, 131B, 1585-1605: critical ice nuclei
        ZNCRIT_REN=5.4E10_JPRB*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB*&
       & (ZSCRITHOMO/(ZSCRITHOMO-1.0_JPRB))**1.5_JPRB/ &
       & (PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB)
        ZNCRIT_REN=ZNCRIT_REN/1.E6_JPRB ! cm**-3

! from Re derive the number concentration - here ice density is 900 kg/m**3 
! Re is in microns, 1e18 factor
        ZCLD=PI(JL,JK)/MAX(PA(JL,JK),ZEPSEC)
        ZCLD=MIN(MAX(ZCLD,0.0_JPRB),ZCLDMAX)
        IF (ZCLD>ZEPSEC) THEN
          ZNICEHOMO=0.75_JPRB*ZRHO(JL,JK)*ZCLD/(RPI*ZRHO_ICE*1.0E-18_JPRB*ZRE_ICE(JL,JK)**3.0_JPRB)
        ENDIF
        ZNICEHOMO = ZNICEHOMO/1.E6_JPRB ! cm**-3

! following Ren and Mackensie, 2005, QJ 131, linearly interpolate to get Ice number
        IF (ZICENUCLEI(JL,JK)<ZNCRIT_REN) THEN
          ZNICE(JL,JK)=ZICENUCLEI(JL,JK)+(1.0_JPRB-ZICENUCLEI(JL,JK)/ZNCRIT_REN)*ZNICEHOMO
        ELSE
          ZNICE(JL,JK)=ZICENUCLEI(JL,JK) 
        ENDIF
        ZDNC2(JL,JK)=ZNICE(JL,JK)

! number is 3/(4*pi*rho_liq*10^6)  [10^6 for N in cm**-3]
        ZRE_ICE(JL,JK)=(0.75_JPRB*ZRHO(JL,JK)*ZCLD/(RPI*ZRHO_ICE*1.E6_JPRB*ZNICE(JL,JK)))**0.333_JPRB
        ZRE_ICE(JL,JK)=ZRE_ICE(JL,JK)*1.E6_JPRB
        ZPDE(JL,JK) = ZRE_ICE(JL,JK)
      ENDIF

    ENDDO
    IF (LLPRINT .AND. LLCLOUD(JK)) THEN
      WRITE(NULOUT,9004) JK,(JL,ZCCN0(JL),ZQLWP(JL,JK),ZRHO(JL,JK),ZDNC2(JL,JK),ZPDE(JL,JK),JL=KIDIA,KFDIA,15)
9004  FORMAT(1X,'aer_diag1 CCN0 LWP RH0 DNC2 ZPDE:',I3,5(2X,I2,5E10.3))
    ENDIF
  ENDDO
!-- ENDIF

!     --------------------------------------------------------------
! Storing the results in PEXTRA field

CALL UPDATE_FIELDS(YDPHY2,2,KIDIA,KFDIA,KLON,KLEV, &
  & PI1 = ZIWC,  PI2 = ZLWC,  PI3 = ZDDE,  PI4 = ZPDE,  PI5 = ZDRE,  PI6 = ZPRE,&
  & PI7 = ZCCND, PI8 = ZCCNP, PI9 = ZDNC1, PI10= ZDNC2, &
  & PO1=PEXTRA(:,:,KIEXT3D+1), PO2=PEXTRA(:,:,KIEXT3D+2), PO3=PEXTRA(:,:,KIEXT3D+3), &
  & PO4=PEXTRA(:,:,KIEXT3D+4), PO5=PEXTRA(:,:,KIEXT3D+5), PO6=PEXTRA(:,:,KIEXT3D+6), &
  & PO7=PEXTRA(:,:,KIEXT3D+7), PO8=PEXTRA(:,:,KIEXT3D+8), PO9=PEXTRA(:,:,KIEXT3D+9), &
  & PO10=PEXTRA(:,:,KIEXT3D+10))

!     --------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_DIAG1',1,ZHOOK_HANDLE)
END SUBROUTINE AER_DIAG1
