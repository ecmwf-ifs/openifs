! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADLSW&
 & (YDDIMV, YDML_PHY_RAD,YDEPHLI,YDPHNC,YDECLD,KIDIA, KFDIA, KLON, KLEV, KMODE, KAER,&
 & PRII0,&
 & PAER , PALBD , PALBP, PAPH , PAP,&
 & PCCNL, PCCNO,&
 & PCCO2, PCLFR , PDP  , PEMIS, PEMIW , PLSM , PMU0, POZON,&
 & PQ   , PQIWP , PQLWP, PQS  ,&
 & PTH  , PT    , PTS  , PNBAS, PNTOP,&
 & PEMIT, PFCT  , PFLT , PFCS , PFLS,&
 & PFRSOD,PSUDU , PUVDF, PPARF, PPARCF, PTINCF,&
 & PSFSWDIR, PSFSWDIF,PFSDNN,PFSDNV,&
 & LDDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST,&
 & PAERINDS)

!**** *RADLSW* - INTERFACE TO ECMWF LW AND SW RADIATION SCHEMES

!     PURPOSE.
!     --------
!           CONTROLS RADIATION COMPUTATIONS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PAERINDS  : (KLON,KLEV)     ; OPTICAL THICKNESS OF THE SULFATE AEROSOL
!                            ; used for indirect effect computation only
! PALBD  : (KLON,NSW)        ; SURF. SW ALBEDO FOR DIFFUSE RADIATION
! PALBP  : (KLON,NSW)        ; SURF. SW ALBEDO FOR PARALLEL RADIATION
! PAPH   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PAP    : (KLON,KLEV)       ; FULL LEVEL PRESSURE
! PCCNL  : (KLON)            ; CCN CONCENTRATION OVER LAND
! PCCNO  : (KLON)            ; CCN CONCENTRATION OVER OCEAN
! PCCO2  :                   ; CONCENTRATION IN CO2 (KG/KG)
! PCLFR  : (KLON,KLEV)       ; CLOUD FRACTIONAL COVER
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PLSM   : (KLON)            ; LAND-SEA MASK
! PMU0   : (KLON)            ; SOLAR ANGLE
! PNBAS  : (KLON)            ; INDEX OF BASE OF CONVECTIVE LAYER
! PNTOP  : (KLON)            ; INDEX OF TOP OF CONVECTIVE LAYER
! POZON  : (KLON,KLEV)       ; OZONE AMOUNT in LAYER (KG/KG*PA)
! PQ     : (KLON,KLEV)       ; SPECIFIC HUMIDITY KG/KG
! PQIWP  : (KLON,KLEV)       ; SOLID  WATER KG/KG
! PQLWP  : (KLON,KLEV)       ; LIQUID WATER KG/KG
! PQS    : (KLON,KLEV)       ; SATURATION WATER VAPOR  KG/KG
! PTH    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
! PT     : (KLON,KLEV)       ; FULL LEVEL TEMPERATURE
! PTS    : (KLON)            ; SURFACE TEMPERATURE
! LDDUST                     ; Dust properties switch
! PPIZA_DST  : (KPROMA,KLEV,NSW); Single scattering albedo of dust 
! PCGA_DST   : (KPROMA,KLEV,NSW); Assymetry factor for dust 
! PTAUREL_DST: (KPROMA,KLEV,NSW); Optical depth of dust relative to at 550nm
!     ==== OUTPUTS ===
! PFCT   : (KLON,KLEV+1)     ; CLEAR-SKY LW NET FLUXES
! PFLT   : (KLON,KLEV+1)     ; TOTAL LW NET FLUXES
! PFCS   : (KLON,KLEV+1)     ; CLEAR-SKY SW NET FLUXES
! PFLS   : (KLON,KLEV+1)     ; TOTAL SW NET FLUXES
! PFRSOD : (KLON)            ; TOTAL-SKY SURFACE SW DOWNWARD FLUX
! PEMIT  : (KLON)            ; SURFACE TOTAL LONGWAVE EMISSIVITY
! PSUDU  : (KLON)            ; SOLAR RADIANCE IN SUN'S DIRECTION
! PPARF  : (KLON)            ; PHOTOSYNTHETICALLY ACTIVE RADIATION
! PUVDF  : (KLON)            ; UV(-B) RADIATION
! PPARCF : (KLON)            ; CLEAR-SKY PHOTOSYNTHETICALLY ACTIVE RADIATION
! PTINCF : (KLON)            ; TOA INCIDENT SOLAR RADIATION 

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHORS.
!     --------
!      J.-J. MORCRETTE         *ECMWF*
!      ORIGINAL : 88-02-04

!     MODIFICATIONS.
!     --------------
!      JJMorcrette : 010112 Sun-Rikus ice particle Diameter
!      JJMorcrette : 010301 cleaning liq/ice cloud optical properties
!      JJMorcrette : 011005 CCN --> Re liquid water clouds
!      JJMorcrette : 011108 Safety checks
!      JJMorcrette : 011108 Safety checks
!      DJSalmond   : 020211 Check before R-To-R
!      JJMorcrette : 020901 PAR & UV
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      JJMorcrette : 050402 New sets of optical properties (NB: inactive)
!      Y.Seity       04-11-18 : add 4 arguments for AROME externalized surface
!      Y.Seity       05-10-10 : add 3 optional arg. for dust SW properties
!      JJMorcrette 20060721 PP of clear-sky PAR and TOA incident solar radiation 
!      M.Janiskova : 22-Nov-2006 Modified call for LW
!      JJMorcrette 20080424 3D input fields for CO2, CH4, N2O, CFC11, 12, 22 and CCL4 (not here)
!        A.Voldoire  2011-02  sulfate indirect effect computation 
!        M.Janiskova : 02-Mar-2012 Initialization of cloud LW properties for RRTM+cleaning
!      P. Lopez     2014-03-24 Retuning of cloud top definition
!                             for computational optimization.
!      P. Lopez     2015-02-11 Added security in cloud top definition
!                              when computational optimization is on.
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOECLD   , ONLY : TECLD
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG       ,RD       ,RTT      ,RPI
USE YOELW    , ONLY : NSIL     ,NTRA     ,NUA      ,TSTAND   ,XP
USE YOERDU   , ONLY : NUAER    ,NTRAER   ,REPLOG   ,REPSC    ,REPSCW   ,DIFF
USE YOETHF   , ONLY : RTICE
USE YOEPHLI  , ONLY : TEPHLI
USE YOERRTWN , ONLY :                     DELWAVE   ,TOTPLNK   
USE YOPHNC   , ONLY : TPHNC

USE YOMLUN   , ONLY : NULOUT
USE YOMCT3   , ONLY : NSTEP

!  Cloud liquid and ice optical properties (N{SW,LW}{LIQ,ICE}OPT
!  Default values set in ../phys_radi/suecrad.F90,
!  modified in namelist &NAERAD

!  liquid water cloud 0: Fouquart    (SW), Smith-Shi   (LW)
!                     1: Slingo      (SW), Savijarvi   (LW)
!                     2: Slingo      (SW), Lindner-Li  (LW)
!                     3: Nielsen     (SW), Smith-Shi   (LW)
!  ice water cloud    0: Ebert-Curry (SW), Smith-Shi   (LW)
!                     1: Ebert-Curry (SW), Ebert-Curry (LW)
!                     2: Fu-Liou'93  (SW), Fu-Liou'93  (LW)
!                     3: Fu'96       (SW), Fu et al'98 (LW)

! Cloud liquid/ice effective radius/diameter (NRADLP, NRADIP)

!  liquid water cloud 0: f(P) 10 to 45
!                     1: 13: ocean; 10: land
!                     2: Martin et al. CCN 50 over ocean, 900 over land
!  ice water cloud    0: 40 microns
!                     1: f(T) 40 to 130 microns
!                     2: f(T) 30 to 60
!                     3: f(T,IWC) Sun'01: 22.5 to 175 microns


IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TECLD)       ,INTENT(IN)    :: YDECLD
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(INOUT):: YDML_PHY_RAD
TYPE(TPHNC)       ,INTENT(INOUT) :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMODE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRII0 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERINDS(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCNL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCNO(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCO2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLFR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZON(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQIWP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLWP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNBAS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV,YDML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV,YDML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV,YDML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOD(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSUDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVDF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARCF(KLON), PTINCF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSFSWDIR(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSFSWDIF(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDNN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDNV(KLON) 

!     -----------------------------------------------------------------

INTEGER(KIND=JPIM) :: IBAS(KLON)     , ITOP(KLON)

REAL(KIND=JPRB) ::&
 & ZALBD(KLON,YDML_PHY_RAD%YRERAD%NSW)    , ZALBP(KLON,YDML_PHY_RAD%YRERAD%NSW)&
 & , ZCG(KLON,YDML_PHY_RAD%YRERAD%NSW,KLEV) , ZOMEGA(KLON,YDML_PHY_RAD%YRERAD%NSW,KLEV)&
 & , ZTAU (KLON,YDML_PHY_RAD%YRERAD%NSW,KLEV)&
 & , ZTAUCLD(KLON,KLEV,16), ZTCLEAR(KLON)  
REAL(KIND=JPRB) ::&
 & ZCLDLD(KLON,KLEV)  , ZCLDLU(KLON,KLEV)&
 & , ZCLDSW(KLON,KLEV)  , ZCLD0(KLON,KLEV)&
 & , ZDT0(KLON)&
 & , ZEMIS(KLON)        , ZEMIW(KLON)&
 & , ZFLUX (KLON,2,KLEV+1)                 , ZFLUC(KLON,2,KLEV+1)&
 & , ZFIWP(KLON)        , ZFLWP(KLON)      , ZFRWP(KLON)&
 & , ZIWC(KLON)         , ZLWC(KLON)&
 !cc            , ZRWC(KLON)
 & , ZMU0(KLON)         , ZOZ(KLON,KLEV)   , ZOZN(KLON,KLEV)&
 & , ZPMB(KLON,KLEV+1)  , ZPSOL(KLON)&
 & , ZTAVE (KLON,KLEV)  , ZTL(KLON,KLEV+1)&
 & , ZVIEW(KLON)  
REAL(KIND=JPRB) ::&
 & ZFCDWN(KLON,KLEV+1), ZFCUP(KLON,KLEV+1)&
 & , ZFSDWN(KLON,KLEV+1), ZFSUP(KLON,KLEV+1)&
 & , ZFSUPN(KLON)       , ZFSUPV(KLON)&
 & , ZFCUPN(KLON)       , ZFCUPV(KLON)&
 & , ZFSDNN(KLON)       , ZFSDNV(KLON)&
 & , ZFCDNN(KLON)       , ZFCDNV(KLON)&
 & , ZDIRFS(KLON,YDML_PHY_RAD%YRERAD%NSW)   , ZDIFFS(KLON,YDML_PHY_RAD%YRERAD%NSW)  
REAL(KIND=JPRB) ::&
 & ZALFICE(KLON)      , ZGAMICE(KLON)     , ZBICE(KLON)   , ZDESR(KLON)&
 & , ZRADIP(KLON)       , ZRADLP(KLON)&
 !cc           , ZRADRD(KLON)
 & , ZRAINT(KLON)       , ZRES(KLON)&
 & , ZTICE(KLON)        , ZEMIT(KLON),  ZBICFU(KLON)&
 & , ZKICFU(KLON)
REAL(KIND=JPRB) :: ZSUDU(KLON)   , ZPARF(KLON)       , ZUVDF(KLON), ZPARCF(KLON)
REAL(KIND=JPRB) :: ZCOND(KLON,KLEV)

REAL(KIND=JPRB) :: ZCO2(KLON,KLEV), ZCH4(KLON,KLEV), ZN2O(KLON,KLEV),&
 & ZNO2(KLON,KLEV), ZC11(KLON,KLEV), ZC12(KLON,KLEV), ZC22(KLON,KLEV), ZCL4(KLON,KLEV)

INTEGER(KIND=JPIM) :: IKL, JK, JKL, JKLP1, JL, JNU, JRTM, JSW, INDLAY

REAL(KIND=JPRB) :: ZASYMX, ZDIFFD, ZGI, ZGL, ZGR, ZIWGKG, ZLWGKG,&
 & ZMSAID, ZMSAIU, ZMSALD, ZMSALU, ZRSAIA, ZRSAID, ZRSAIE, ZRSAIF, ZRSAIG, ZRSALD,&
 & ZMULTI, ZMULTL, ZOI   , ZOL,&
 & ZOMGMX, ZOR, ZRMUZ, ZRWGKG, ZTAUMX, ZTEMPC,&
 & ZTOI, ZTOL, ZTOR, ZZFIWP, ZZFLWP, ZDPOG, ZPODT  

REAL(KIND=JPRB) :: ZALND, ZASEA, ZD, ZDEN, ZNTOT, ZNUM, ZRATIO, Z1RADI,&
 & Z1RADL, ZBETAI, ZOMGI, ZOMGP, ZFDEL, ZTCELS, ZFSR, ZAIWC,&
 & ZBIWC, ZTBLAY, ZADDPLK, ZPLANCK, ZEXTCF, Z1MOMG,&
 & ZDEFRE, ZREFDE, ZVI , ZMABSD 

!REAL(KIND=JPRB) :: ZAVDP(KLON), ZAVTO(KLON), ZSQTO(KLON)
REAL(KIND=JPRB) :: ZAVTO(KLON), ZSQTO(KLON)
REAL(KIND=JPRB) :: ZSQUAR(KLON,KLEV), ZVARIA(KLON,KLEV)
INTEGER(KIND=JPIM) :: IKI, JKI, JEXPLR, JXPLDN
INTEGER(KIND=JPIM) :: IUAER, ITOPC, ILEV(KLON), ILEVR

REAL(KIND=JPRB) :: ZBICFU_INIT
INTEGER(KIND=JPIM) :: ICE_WATER_CLOUD_TOT_EMIS
REAL(KIND=JPRB) :: ZDZHT(KLON,KLEV), ZMAER(KLON,KLEV)
REAL(KIND=JPRB) :: ZCCDNC(KLON,KLEV), ZCDNC(KLON,KLEV), ZRAVRG(KLON,KLEV)
REAL(KIND=JPRB) :: ZTAUTAER, ZALPHAER, ZCOEFRFH


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "lw.intfb.h"
#include "rrtm_rrtm_140gp.intfb.h"
#include "sw.intfb.h"
#include "abor1.intfb.h"

!     -----------------------------------------------------------------

!*         1.     SET-UP INPUT QUANTITIES FOR RADIATION
!                 -------------------------------------

IF (LHOOK) CALL DR_HOOK('RADLSW',0,ZHOOK_HANDLE)
ASSOCIATE(LWLCLHR=>YDML_PHY_RAD%YRELWRAD%LWLCLHR, NLEVLWC=>YDML_PHY_RAD%YRELWRAD%NLEVLWC, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, &
 & LCCNL=>YDML_PHY_RAD%YRERAD%LCCNL, LCCNO=>YDML_PHY_RAD%YRERAD%LCCNO, LDIFFC=>YDML_PHY_RAD%YRERAD%LDIFFC, &
 & LEDBUG=>YDML_PHY_RAD%YRERAD%LEDBUG, LNEWAER=>YDML_PHY_RAD%YRERAD%LNEWAER, LRRTM=>YDML_PHY_RAD%YRERAD%LRRTM, &
 & NICEOPT=>YDML_PHY_RAD%YRERAD%NICEOPT, NINHOM=>YDML_PHY_RAD%YRERAD%NINHOM, NLAYINH=>YDML_PHY_RAD%YRERAD%NLAYINH, &
 & NLIQOPT=>YDML_PHY_RAD%YRERAD%NLIQOPT, NRADIP=>YDML_PHY_RAD%YRERAD%NRADIP, NRADLP=>YDML_PHY_RAD%YRERAD%NRADLP, &
 & NSW=>YDML_PHY_RAD%YRERAD%NSW, RCCNLND=>YDML_PHY_RAD%YRERAD%RCCNLND, RCCNSEA=>YDML_PHY_RAD%YRERAD%RCCNSEA, &
 & RLWINHF=>YDML_PHY_RAD%YRERAD%RLWINHF, RRE2DE=>YDML_PHY_RAD%YRERAD%RRE2DE, RSWINHF=>YDML_PHY_RAD%YRERAD%RSWINHF, &
 & NSWICEOPT=>YDML_PHY_RAD%YRERAD%NSWICEOPT, NSWLIQOPT=>YDML_PHY_RAD%YRERAD%NSWLIQOPT, &
 & NLWICEOPT=>YDML_PHY_RAD%YRERAD%NLWICEOPT, NLWLIQOPT=>YDML_PHY_RAD%YRERAD%NLWLIQOPT, &
 & RCCL4=>YDML_PHY_RAD%YRERDI%RCCL4, RCFC11=>YDML_PHY_RAD%YRERDI%RCFC11, RCFC12=>YDML_PHY_RAD%YRERDI%RCFC12, &
 & RCFC22=>YDML_PHY_RAD%YRERDI%RCFC22, RCH4=>YDML_PHY_RAD%YRERDI%RCH4, RN2O=>YDML_PHY_RAD%YRERDI%RN2O, &
 & LH2OCO2=>YDPHNC%LH2OCO2, &
 & RYFWCA=>YDML_PHY_RAD%YRESWRT%RYFWCA ,RYFWCB=>YDML_PHY_RAD%YRESWRT%RYFWCB ,RYFWCC=>YDML_PHY_RAD%YRESWRT%RYFWCC ,&
 & RYFWCD=>YDML_PHY_RAD%YRESWRT%RYFWCD ,RYFWCE=>YDML_PHY_RAD%YRESWRT%RYFWCE ,RYFWCF=>YDML_PHY_RAD%YRESWRT%RYFWCF ,& 
 & REBCUA=>YDML_PHY_RAD%YRESWRT%REBCUA ,REBCUB=>YDML_PHY_RAD%YRESWRT%REBCUB ,REBCUC=>YDML_PHY_RAD%YRESWRT%REBCUC ,&
 & REBCUD=>YDML_PHY_RAD%YRESWRT%REBCUD ,REBCUE=>YDML_PHY_RAD%YRESWRT%REBCUE ,REBCUF=>YDML_PHY_RAD%YRESWRT%REBCUF ,&
 & REBCUI=>YDML_PHY_RAD%YRESWRT%REBCUI ,REBCUJ=>YDML_PHY_RAD%YRESWRT%REBCUJ ,REBCUG=>YDML_PHY_RAD%YRESWRT%REBCUG ,&
 & REBCUH=>YDML_PHY_RAD%YRESWRT%REBCUH ,RHSAVI=>YDML_PHY_RAD%YRESWRT%RHSAVI ,RFULIO=>YDML_PHY_RAD%YRESWRT%RFULIO ,&
 & RFLAA0=>YDML_PHY_RAD%YRESWRT%RFLAA0 ,RFLAA1=>YDML_PHY_RAD%YRESWRT%RFLAA1 ,RFLBB0=>YDML_PHY_RAD%YRESWRT%RFLBB0 ,&
 & RFLBB1=>YDML_PHY_RAD%YRESWRT%RFLBB1 ,RFLBB2=>YDML_PHY_RAD%YRESWRT%RFLBB2 ,RFLBB3=>YDML_PHY_RAD%YRESWRT%RFLBB3 ,&
 & RFLCC0=>YDML_PHY_RAD%YRESWRT%RFLCC0 ,RFLCC1=>YDML_PHY_RAD%YRESWRT%RFLCC1 ,RFLCC2=>YDML_PHY_RAD%YRESWRT%RFLCC2 ,&
 & RFLCC3=>YDML_PHY_RAD%YRESWRT%RFLCC3 ,RFLDD0=>YDML_PHY_RAD%YRESWRT%RFLDD0 ,RFLDD1=>YDML_PHY_RAD%YRESWRT%RFLDD1 ,&
 & RFLDD2=>YDML_PHY_RAD%YRESWRT%RFLDD2 ,RFLDD3=>YDML_PHY_RAD%YRESWRT%RFLDD3 ,RFUETA=>YDML_PHY_RAD%YRESWRT%RFUETA ,&
 & RFUETB=>YDML_PHY_RAD%YRESWRT%RFUETB ,RFUETC=>YDML_PHY_RAD%YRESWRT%RFUETC ,RASWCA=>YDML_PHY_RAD%YRESWRT%RASWCA ,&
 & RASWCB=>YDML_PHY_RAD%YRESWRT%RASWCB ,RASWCC=>YDML_PHY_RAD%YRESWRT%RASWCC ,RASWCD=>YDML_PHY_RAD%YRESWRT%RASWCD ,&
 & RASWCE=>YDML_PHY_RAD%YRESWRT%RASWCE ,RASWCF=>YDML_PHY_RAD%YRESWRT%RASWCF ,RKPNWCA=>YDML_PHY_RAD%YRESWRT%RKPNWCA,&
 & RKPNWCB=>YDML_PHY_RAD%YRESWRT%RKPNWCB,RKPNWCC=>YDML_PHY_RAD%YRESWRT%RKPNWCC,RKPNWCD=>YDML_PHY_RAD%YRESWRT%RKPNWCD,&
 & RKPNWCE=>YDML_PHY_RAD%YRESWRT%RKPNWCE,RKPNWCF=>YDML_PHY_RAD%YRESWRT%RKPNWCF,RKPNWCG=>YDML_PHY_RAD%YRESWRT%RKPNWCG,&
 & RKPNWCH=>YDML_PHY_RAD%YRESWRT%RKPNWCH,RFUAA0=>YDML_PHY_RAD%YRESWRT%RFUAA0 ,RFUAA1=>YDML_PHY_RAD%YRESWRT%RFUAA1,&
 & RFUBB0=>YDML_PHY_RAD%YRESWRT%RFUBB0 ,RFUBB1=>YDML_PHY_RAD%YRESWRT%RFUBB1 ,RFUBB2=>YDML_PHY_RAD%YRESWRT%RFUBB2 ,&
 & RFUBB3=>YDML_PHY_RAD%YRESWRT%RFUBB3 ,RFUCC0=>YDML_PHY_RAD%YRESWRT%RFUCC0 ,RFUCC1=>YDML_PHY_RAD%YRESWRT%RFUCC1 ,&
 & RFUCC2=>YDML_PHY_RAD%YRESWRT%RFUCC2 ,RFUCC3=>YDML_PHY_RAD%YRESWRT%RFUCC3 ,&
 & RLILIA=>YDML_PHY_RAD%YRESWRT%RLILIA ,RLILIB=>YDML_PHY_RAD%YRESWRT%RLILIB  )

ZREFDE = RRE2DE
ZDEFRE = 1.0_JPRB / ZREFDE

DO JL = KIDIA,KFDIA
  ZFCUP(JL,KLEV+1) = 0.0_JPRB
  ZFCDWN(JL,KLEV+1) = REPLOG
  ZFSUP(JL,KLEV+1) = 0.0_JPRB
  ZFSDWN(JL,KLEV+1) = REPLOG
  ZFLUX(JL,1,KLEV+1) = 0.0_JPRB
  ZFLUX(JL,2,KLEV+1) = 0.0_JPRB
  ZFLUC(JL,1,KLEV+1) = 0.0_JPRB
  ZFLUC(JL,2,KLEV+1) = 0.0_JPRB
  ZFSDNN(JL) = 0.0_JPRB
  ZFSDNV(JL) = 0.0_JPRB
  ZFCDNN(JL) = 0.0_JPRB
  ZFCDNV(JL) = 0.0_JPRB
  ZFSUPN(JL) = 0.0_JPRB
  ZFSUPV(JL) = 0.0_JPRB
  ZFCUPN(JL) = 0.0_JPRB
  ZFCUPV(JL) = 0.0_JPRB
  ZPSOL(JL) = PAPH(JL,KLEV+1)
  ZPMB(JL,1) = ZPSOL(JL) / 100.0_JPRB
  ZDT0(JL) = PTS(JL) - PTH(JL,KLEV+1)
  PSUDU(JL) = 0.0_JPRB
  PPARF(JL) = 0.0_JPRB
  PPARCF(JL)= 0.0_JPRB
  PUVDF(JL) = 0.0_JPRB
  PSFSWDIR(JL,:)=0.0_JPRB
  PSFSWDIF(JL,:)=0.0_JPRB
  IBAS(JL) = INT ( 0.01_JPRB + PNBAS(JL) )
  ITOP(JL) = INT ( 0.01_JPRB + PNTOP(JL) )
ENDDO

!*         1.1    INITIALIZE VARIOUS FIELDS
!                 -------------------------

DO JSW=1,NSW
  DO JL = KIDIA,KFDIA
    ZALBD(JL,JSW)=PALBD(JL,JSW)
    ZALBP(JL,JSW)=PALBP(JL,JSW)
  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZEMIS(JL)  =PEMIS(JL)
  ZEMIW(JL)  =PEMIW(JL)
  ZMU0(JL)   =PMU0(JL)
ENDDO

DO JK = 1 , KLEV
  JKL = KLEV+ 1 - JK
  DO JL = KIDIA,KFDIA
    ZPMB(JL,JK+1)=PAPH(JL,JKL)/100.0_JPRB

!-- ZOZ in cm.atm for SW scheme    
    ZOZ(JL,JK)   = POZON(JL,JKL) * 46.6968_JPRB / RG

    ZCLD0(JL,JK) = 0.0_JPRB
    ZFCUP(JL,JK) = 0.0_JPRB
    ZFCDWN(JL,JK) = 0.0_JPRB
    ZFSUP(JL,JK) = 0.0_JPRB
    ZFSDWN(JL,JK) = 0.0_JPRB
    ZFLUX(JL,1,JK) = 0.0_JPRB
    ZFLUX(JL,2,JK) = 0.0_JPRB
    ZFLUC(JL,1,JK) = 0.0_JPRB
    ZFLUC(JL,2,JK) = 0.0_JPRB
  ENDDO
ENDDO

DO JK = 1 , KLEV
  DO JL = KIDIA,KFDIA
    ZTAUCLD(JL,JK,:) = 0.0_JPRB
  ENDDO
ENDDO

DO JK=1,KLEV
  JKL=KLEV+1-JK
  JKLP1=JKL+1
  DO JL=KIDIA,KFDIA
    ZTL(JL,JK)=PTH(JL,JKLP1)
    ZTAVE(JL,JK)=PT(JL,JKL)
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  ZTL(JL,KLEV+1)= PTH(JL,1)
  ZPMB(JL,KLEV+1) = PAPH(JL,1)/100.0_JPRB
ENDDO
!***

!     ------------------------------------------------------------------

!*         2.     CLOUD AND AEROSOL PARAMETERS
!                 ----------------------------

DO JK = 1 , KLEV
  IKL = KLEV + 1 - JK

!          2.1    INITIALIZE OPTICAL PROPERTIES TO CLEAR SKY VALUES
!                 -------------------------------------------------

  DO JSW = 1,NSW
    DO JL = KIDIA,KFDIA
      ZTAU(JL,JSW,JK)  = 0.0_JPRB
      ZOMEGA(JL,JSW,JK)= 1.0_JPRB
      ZCG(JL,JSW,JK)   = 0.0_JPRB
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZCLDSW(JL,JK)  = 0.0_JPRB
    ZCLDLD(JL,JK)  = 0.0_JPRB
    ZCLDLU(JL,JK)  = 0.0_JPRB
  ENDDO

!          2.2    CLOUD ICE AND LIQUID CONTENT AND PATH
!                 -------------------------------------

  DO JL = KIDIA,KFDIA

! --- LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
    IF (PCLFR(JL,IKL) > REPSC ) THEN
      ZLWGKG=MAX(PQLWP(JL,IKL)*1000.0_JPRB,0.0_JPRB)
      ZIWGKG=MAX(PQIWP(JL,IKL)*1000.0_JPRB,0.0_JPRB)
      ZLWGKG=ZLWGKG/PCLFR(JL,IKL)
      ZIWGKG=ZIWGKG/PCLFR(JL,IKL)
    ELSE
      ZLWGKG=0.0_JPRB
      ZIWGKG=0.0_JPRB
    ENDIF
    ZRWGKG=0.0_JPRB
    ZRAINT(JL)=0.0_JPRB

! --- RAIN LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
!    IF (PRAINT(JL,IKL) >= REPSCW) THEN
!      ZRWGKG=MAX(PQRAIN(JL,IKL)*1000., 0.0)
!      ZRAINT(JL)=PRAINT(JL,IKL)*3600.*1000.
!- no radiative effect of rain (for the moment)
!      ZRWGKG=0.
!      ZRAINT(JL)=0.
! ===========================================================

! Modifications Martin et al.
!    ELSE
!    ENDIF
    ZDPOG=PDP(JL,IKL)/RG
    ZFLWP(JL)= ZLWGKG*ZDPOG
    ZFIWP(JL)= ZIWGKG*ZDPOG
    ZFRWP(JL)= ZRWGKG*ZDPOG
    ZPODT=PAP(JL,IKL)/(RD*PT(JL,IKL))
    ZLWC(JL)=ZLWGKG*ZPODT
    ZIWC(JL)=ZIWGKG*ZPODT
!    ZRWC(JL)=ZRWGKG*ZPODT

  ENDDO

  DO JL = KIDIA,KFDIA
! --- EFFECTIVE RADIUS FOR WATER, ICE AND RAIN PARTICLES

    IF (NRADLP == 9) THEN
! --- fixed re for liquid drops (for testing purposes)
      ZRADLP(JL)=10.0_JPRB

    ELSEIF (NRADLP == 0) THEN
!-- very old parametrization as f(pressure) ERA-15
      ZRADLP(JL)=10.0_JPRB + (100000.0_JPRB-PAP(JL,IKL))*3.5_JPRB

    ELSEIF (NRADLP == 1) THEN
! simple distinction between land (10) and ocean (13) Zhang and Rossow
      IF (PLSM(JL) < 0.5_JPRB) THEN
        ZRADLP(JL)=13.0_JPRB
      ELSE
        ZRADLP(JL)=10.0_JPRB
      ENDIF
      
    ELSEIF (NRADLP == 2) THEN
!--  based on Martin et al., 1994, JAS
      IF (PLSM(JL) < 0.5_JPRB) THEN
        IF (LCCNO) THEN
!          ZASEA=50.0_JPRB
          ZASEA=PCCNO(JL)
        ELSE  
          ZASEA=RCCNSEA
        ENDIF  
        ZD=0.33_JPRB
        ZNTOT=-1.15E-03_JPRB*ZASEA*ZASEA+0.963_JPRB*ZASEA+5.30_JPRB
      ELSE
        IF (LCCNL) THEN 
!          ZALND=900.0_JPRB
          ZALND=PCCNL(JL)
        ELSE  
          ZALND=RCCNLND
        ENDIF  
        ZD=0.43_JPRB
        ZNTOT=-2.10E-04_JPRB*ZALND*ZALND+0.568_JPRB*ZALND-27.9_JPRB
      ENDIF
      ZNUM=3.0_JPRB*ZLWC(JL)*(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
      ZDEN=4.0_JPRB*RPI*ZNTOT*(1.0_JPRB+ZD*ZD)**3
      IF((ZNUM/ZDEN) > REPLOG)THEN
        ZRADLP(JL)=100.0_JPRB*EXP(0.333_JPRB*LOG(ZNUM/ZDEN))
        ZRADLP(JL)=MAX(ZRADLP(JL), 4.0_JPRB)
        ZRADLP(JL)=MIN(ZRADLP(JL),16.0_JPRB)
      ELSE
        ZRADLP(JL)=4.0_JPRB
      ENDIF
! Parametrisation pour l'effet indirect des aerosols sulfates      
    ELSEIF(NRADLP == 3) THEN
      IF (.NOT. LNEWAER) THEN
       CALL ABOR1('RADLSW: CASE NOT CODED NRADLP=3/LNEWAER=.T.')
      ELSE
! --- EFFECTIVE RADIUS FOR WATER AND ICE PARTICLES

        ZTAUTAER=PAERINDS(JL,IKL)

!+    for indirect effect of aerosols, HRM, Nov.24, 1998.

        ZALPHAER=5.0_JPRB
        ZCOEFRFH=1.7_JPRB
        ZDZHT(JL,IKL)=RD*PT(JL,IKL)*PDP(JL,IKL)/(PAP(JL,IKL)*RG)
        ZMAER(JL,IKL)=1.0E6_JPRB*ZTAUTAER/(ZALPHAER*ZCOEFRFH*ZDZHT(JL,IKL))
!
! .. This is to avoid a problem when taking the log as CDNC is
! computed. Has no influence on CDNC even if ZMAER<1e-10, as 
! CDNC is taken to be MAX(CDNC,20) afterwards.
!
        ZMAER(JL,IKL)=MAX(ZMAER(JL,IKL),1.0E-10_JPRB)
!
! .. Modification of Boucher(1994)'s formula giving the number of 
! condensation nuclei. The coefficients were obtained from a calibration
! on POLDER measurements (J. Quaas, personal communication).
!
        ZCCDNC(JL,IKL)=10.0_JPRB**(1.7_JPRB+0.2_JPRB*LOG10(ZMAER(JL,IKL)))
!
! .. Threshold on the number of droplets (O. Boucher, J. Quaas, personal
! communication)
!
        ZCCDNC(JL,IKL)=MAX(ZCCDNC(JL,IKL),20.0_JPRB)
        ZCDNC(JL,IKL)=4.0_JPRB*3.14159_JPRB*ZCCDNC(JL,IKL)*1.0E6_JPRB
        ZRAVRG(JL,IKL)=3.0_JPRB*ZFLWP(JL)/ZDZHT(JL,IKL)/ZCDNC(JL,IKL)
!
! .. Calibrate parameterization due to the low modeled liquid water
! content in the model (O. Boucher, personal communication)
!    We also impose that the maximum radius of aerosol droplets 
! is 45 micrometers, a reminder of previous versions, where the 
! radii of these droplets used to be:
!    10 + 0.035*(1000.-P), P=pressure in hPa.
!
        ZRADLP(JL)=1.1_JPRB*ZRAVRG(JL,IKL)**(1.0_JPRB/3.0_JPRB)*1.0E4_JPRB
! Keep min as in previous case
! Keep max as in previous case
        ZRADLP(JL)=MAX(ZRADLP(JL),4.0_JPRB)
        ZRADLP(JL)=MIN(ZRADLP(JL),16.0_JPRB)

      ENDIF
    ENDIF


! ===========================================================
! ___________________________________________________________

! rain drop from          : unused as ZRAINT is 0.
!    ZRADRD(JL)=500.0_JPRB*ZRAINT(JL)**0.22_JPRB
!    IF (ZFLWP(JL).GT.0.) THEN
!      ZRADRD(JL)=ZRADLP(JL)+ZRADRD(JL)
!    ENDIF   

  ENDDO

  DO JL = KIDIA,KFDIA

! diagnosing the ice particle effective radius/diameter

!- ice particle effective radius =f(T) from Liou and Ou (1994)
 
    IF (PT(JL,IKL) < RTICE) THEN
      ZTEMPC=PT(JL,IKL)-RTT
    ELSE
      ZTEMPC=RTICE-RTT
    ENDIF
    ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
      & 0.0012_JPRB))    

    IF (NRADIP == 0) THEN
!-- fixed 50 micron effective radius
      ZRADIP(JL)= 50.0_JPRB
      ZDESR(JL) = ZDEFRE * ZRADIP(JL)
      
    ELSEIF (NRADIP == 1) THEN 

!-- old formulation based on Liou & Ou (1994) temperature (40-130microns)    
      ZRADIP(JL)=MAX(ZRADIP(JL),40.0_JPRB)
      ZDESR(JL) = ZDEFRE * ZRADIP(JL)
      
    ELSEIF (NRADIP == 2) THEN  
!-- formulation following Jakob, Klein modifications to ice content    
      ZRADIP(JL)=MAX(ZRADIP(JL),30.0_JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),60.0_JPRB)
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
 
    ELSEIF (NRADIP == 3  ) THEN
 
!- ice particle effective radius =f(T,IWC) from Sun and Rikus (1999)
! revised by Sun (2001)
      IF (ZIWC(JL) > 0.0_JPRB ) THEN
        ZTEMPC = PT(JL,IKL)-83.15_JPRB
        ZTCELS = PT(JL,IKL)-RTT
        ZFSR = 1.2351_JPRB +0.0105_JPRB * ZTCELS
! Sun, 2001 (corrected from Sun & Rikus, 1999)
        ZAIWC = 45.8966_JPRB * ZIWC(JL)**0.2214_JPRB
        ZBIWC = 0.7957_JPRB * ZIWC(JL)**0.2535_JPRB
        ZDESR(JL) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)
        ZDESR(JL) = MIN ( MAX( ZDESR(JL), 30.0_JPRB), 155.0_JPRB) ! new
        ! ZDESR(JL) = MIN ( MAX( ZDESR(JL), 45.0_JPRB), 350.0_JPRB) ! old
        ZRADIP(JL)= ZREFDE * ZDESR(JL)
        ZRADIP(JL)= ZRADIP(JL) + 15._JPRB
      ELSE
!        ZDESR(JL) = 92.5_JPRB
        ZDESR(JL) = 80.0_JPRB
        ZRADIP(JL)= ZREFDE * ZDESR(JL)
      ENDIF  
    ENDIF  
    
  ENDDO

!          2.3    CLOUD SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------------

!   -------------------------
! --+ SW OPTICAL PARAMETERS +  Water clouds after Fouquart (1987)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

  DO JSW=1,NSW
    DO JL = KIDIA,KFDIA
      ZTOL=0.0_JPRB
      ZGL =0.0_JPRB
      ZOL =0.0_JPRB
      ZTOI=0.0_JPRB
      ZGI =0.0_JPRB
      ZOI =0.0_JPRB
      ZTOR=0.0_JPRB
      ZGR =0.0_JPRB
      ZOR =0.0_JPRB
      IF (ZFLWP(JL)+ZFIWP(JL)+ZFRWP(JL) > 2.0_JPRB * REPSCW ) THEN
        IF (ZFLWP(JL) >= REPSCW ) THEN
          IF (NSWLIQOPT == 0 ) THEN
!-- SW: Fouquart, 1991
            ZTOL = ZFLWP(JL)*(RYFWCA(JSW)+RYFWCB(JSW)/ZRADLP(JL))
            ZGL  = RYFWCF(JSW)
!            ZOL  = RYFWCC(JSW)-RYFWCD(JSW)*EXP(-RYFWCE(JSW)*ZTOL)
!-- NB: RSWINHF is there simply for making the CY29R2 branch bit compatible with 
! the previous. Should be cleaned when RRTM_SW becomes active
            ZOL  = RYFWCC(JSW)-RYFWCD(JSW)*EXP(-RYFWCE(JSW)*ZTOL*RSWINHF)
          ELSE IF (NSWLIQOPT == 3) THEN
!-- SW: Nielsen, 2013
            ZTOL = ZFLWP(JL)*RKPNWCA(JSW)*ZRADLP(JL)**(-RKPNWCB(JSW))
            ZGL  = RKPNWCE(JSW) + RKPNWCF(JSW)*ZRADLP(JL) -&
                 & RKPNWCG(JSW)*EXP(-RKPNWCH(JSW)*ZRADLP(JL))
            ZOL  = RKPNWCC(JSW) - RKPNWCD(JSW)*ZRADLP(JL)
          ELSE    
!-- SW: Slingo, 1989
            ZTOL = ZFLWP(JL)*(RASWCA(JSW)+RASWCB(JSW)/ZRADLP(JL))
            ZGL  = RASWCE(JSW)+RASWCF(JSW)*ZRADLP(JL)
            ZOL  = 1. - RASWCC(JSW)-RASWCD(JSW)*ZRADLP(JL)
          ENDIF 
        ENDIF

        IF (ZFIWP(JL) >= REPSCW ) THEN
          IF (NSWICEOPT <= 1) THEN
!-- SW: Ebert-Curry          
            ZTOI = ZFIWP(JL)*(REBCUA(JSW)+REBCUB(JSW)/ZRADIP(JL))
            ZGI  = REBCUE(JSW)+REBCUF(JSW)*ZRADIP(JL)
            ZOI  = 1.0_JPRB - REBCUC(JSW)-REBCUD(JSW)*ZRADIP(JL)
            
          ELSEIF (NSWICEOPT == 2) THEN  
!-- SW: Fu-Liou 1993
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZBETAI = RFLAA0(JSW)+Z1RADI* RFLAA1(JSW)
            ZTOI = ZFIWP(JL) * ZBETAI
            ZOMGI= RFLBB0(JSW)+ZRADIP(JL)*(RFLBB1(JSW) + ZRADIP(JL)&
             & *(RFLBB2(JSW)+ZRADIP(JL)* RFLBB3(JSW) ))              
            ZOI  = 1.0_JPRB - ZOMGI
            ZOMGP= RFLCC0(JSW)+ZRADIP(JL)*(RFLCC1(JSW) + ZRADIP(JL)&
             & *(RFLCC2(JSW)+ZRADIP(JL)* RFLCC3(JSW) ))   
            ZFDEL= RFLDD0(JSW)+ZRADIP(JL)*(RFLDD1(JSW) + ZRADIP(JL)&
             & *(RFLDD2(JSW)+ZRADIP(JL)* RFLDD3(JSW) ))   
            ZGI  = ((1.0_JPRB -ZFDEL)*ZOMGP + ZFDEL*3.0_JPRB) / 3.0_JPRB
            
          ELSEIF (NSWICEOPT == 3) THEN  
!-- SW: Fu 1996
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZBETAI = RFUAA0(JSW)+Z1RADI* RFUAA1(JSW)
            ZTOI = ZFIWP(JL) * ZBETAI
            ZOMGI= RFUBB0(JSW)+ZDESR(JL)*(RFUBB1(JSW) + ZDESR(JL)&
             &   *(RFUBB2(JSW)+ZDESR(JL)* RFUBB3(JSW) ))            
            ZOI  = 1.0_JPRB - ZOMGI
            ZGI  = RFUCC0(JSW)+ZDESR(JL)*(RFUCC1(JSW) + ZDESR(JL)&
             &   *(RFUCC2(JSW)+ZDESR(JL)* RFUCC3(JSW) )) 
            ZGI  = MIN(1.0_JPRB, ZGI)
           ELSE
              WRITE(NULOUT,*) 'UNDEFINED NSWLIQOPT =',NSWLIQOPT,' NSWICEOPT =',NSWICEOPT
           ENDIF
   
        ENDIF

!        IF (ZFRWP(JL) >= REPSCW ) THEN
!          ZTOR= ZFRWP(JL)*0.003_JPRB * ZRAINT(JL)**(-0.22_JPRB)         
!          ZOR = 1.0_JPRB - RROMA(JSW)*ZRAINT(JL)**RROMB(JSW)
!          ZGR = RRASY(JSW)
!        ENDIF   

!  - MIX of WATER and ICE CLOUDS
        ZTAUMX= ZTOL + ZTOI + ZTOR
        ZOMGMX= ZTOL*ZOL + ZTOI*ZOI + ZTOR*ZOR
        ZASYMX= ZTOL*ZOL*ZGL + ZTOI*ZOI*ZGI + ZTOR*ZOR*ZGR

        ZASYMX= ZASYMX/ZOMGMX
        ZOMGMX= ZOMGMX/ZTAUMX

! --- SW FINAL CLOUD OPTICAL PARAMETERS

        ZCLDSW(JL,JK)  = PCLFR(JL,IKL)
        ZTAU(JL,JSW,JK)  = ZTAUMX
        ZOMEGA(JL,JSW,JK)= ZOMGMX
        ZCG(JL,JSW,JK)   = ZASYMX
      ENDIF
    ENDDO
  ENDDO

!          2.4    CLOUD LONGWAVE OPTICAL PROPERTIES FOR EC-OPE
!                 --------------------------------------------

!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Smith and Shi (1992)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

  IF (.NOT.LRRTM) THEN

    IF (NLWICEOPT == 1) THEN
      ZBICFU_INIT = 1.0_JPRB
    ELSE
      ZBICFU_INIT = 0.0_JPRB
    ENDIF
    DO JL = KIDIA,KFDIA
      ZALFICE(JL)=0.0_JPRB
      ZGAMICE(JL)=0.0_JPRB
      ZBICE(JL)=0.0_JPRB
      ZTICE(JL)=(PT(JL,IKL)-TSTAND)/TSTAND
      ZBICFU(JL)=ZBICFU_INIT
      ZKICFU(JL)=0.0_JPRB
    ENDDO
    
    DO JNU= 1,NSIL
      DO JL = KIDIA,KFDIA
        ZRES(JL)  = XP(1,JNU)+ZTICE(JL)*(XP(2,JNU)+ZTICE(JL)*(XP(3,&
         & JNU)&
         & +ZTICE(JL)*(XP(4,JNU)+ZTICE(JL)*(XP(5,JNU)+ZTICE(JL)*(XP(6,&
         & JNU)&
         & )))))  
        ZBICE(JL) = ZBICE(JL) + ZRES(JL)
        ZGAMICE(JL) = ZGAMICE(JL) + REBCUI(JNU)*ZRES(JL)
        ZALFICE(JL) = ZALFICE(JL) + REBCUJ(JNU)*ZRES(JL)
      ENDDO
    ENDDO
    
!-- Fu et al. (1998) with M'91 LW scheme    
    IF (NLWICEOPT == 2 .OR. NLWICEOPT == 3) THEN
      DO JRTM=1,16
        DO JL=KIDIA,KFDIA
          IF (PT(JL,IKL) < 160.0_JPRB) THEN
            INDLAY=1
            ZTBLAY =PT(JL,IKL)-160.0_JPRB
          ELSEIF (PT(JL,IKL) < 339.0_JPRB ) THEN
            INDLAY=PT(JL,IKL)-159.0_JPRB
            INDLAY=MAX(INDLAY,1)
            ZTBLAY =PT(JL,IKL)-INT(PT(JL,IKL))
          ELSE 
            INDLAY=180
            ZTBLAY =PT(JL,IKL)-339.0_JPRB
          ENDIF
          ZADDPLK = TOTPLNK(INDLAY+1,JRTM)-TOTPLNK(INDLAY,JRTM)
          ZPLANCK = DELWAVE(JRTM) * (TOTPLNK(INDLAY,JRTM) + ZTBLAY*ZADDPLK)
          ZBICFU(JL) = ZBICFU(JL) + ZPLANCK
        
          IF (ZIWC(JL) > 0.0_JPRB ) THEN
            ZRATIO =  1.0_JPRB / ZDESR(JL) 
            IF (NLWICEOPT == 2) THEN
! ice cloud spectral emissivity a la Fu & Liou (1993)
              ZMABSD = RFULIO(JRTM,1) + ZRATIO&
               & *(RFULIO(JRTM,2) + ZRATIO*RFULIO(JRTM,3))  
          
! ice cloud spectral emissivity a la Fu et al (1998)
            ELSEIF (NLWICEOPT == 3) THEN 
              ZMABSD = RFUETA(JRTM,1) + ZRATIO&
               & *(RFUETA(JRTM,2) + ZRATIO*RFUETA(JRTM,3))  
            ENDIF
            ZKICFU(JL) = ZKICFU(JL)+ ZMABSD*ZPLANCK
          ENDIF  
        ENDDO
      ENDDO
    ENDIF
    
    DO JL = KIDIA,KFDIA
      ZGAMICE(JL) = ZGAMICE(JL) / ZBICE(JL)
      ZALFICE(JL) = ZALFICE(JL) / ZBICE(JL)
      ZKICFU(JL)  = ZKICFU(JL) / ZBICFU(JL)
      
      IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN

        IF (NLWLIQOPT == 0 .OR. NLWLIQOPT >= 3 ) THEN
! water cloud emissivity a la Smith & Shi (1992)
          ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
          ZMSALD= 0.158_JPRB*ZMULTL
          ZMSALU= 0.130_JPRB*ZMULTL
          
        ELSE
! water cloud emissivity a la Savijarvi (1997)
          ZMSALU= 0.2441_JPRB-0.0105_JPRB*ZRADLP(JL)
          ZMSALD= 1.2154_JPRB*ZMSALU
          
        ENDIF  
          
        IF (NLWICEOPT == 0) THEN          
! ice cloud emissivity a la Smith & Shi (1992)
          ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
          ZMSAID= 0.113_JPRB*ZMULTI
          ZMSAIU= 0.093_JPRB*ZMULTI

        ELSEIF (NLWICEOPT == 1) THEN
! ice cloud emissivity a la Ebert & Curry (1992)
          ZMSAID= 1.66_JPRB*(ZALFICE(JL)+ZGAMICE(JL)/ZRADIP(JL))
          ZMSAIU= ZMSAID
          
        ELSEIF (NLWICEOPT == 2 .OR. NLWICEOPT == 3) THEN  
! ice cloud emissivity a la Fu & Liou (1993) or Fu et al. (1998)
          ZMSAID= 1.66_JPRB*ZKICFU(JL)
          ZMSAIU= ZMSAID          
        ENDIF 
       
        IF (NINHOM == 1) THEN
          ZZFLWP= ZFLWP(JL) * RLWINHF
          ZZFIWP= ZFIWP(JL) * RLWINHF
        ELSE
          ZZFLWP= ZFLWP(JL)
          ZZFIWP= ZFIWP(JL)
        ENDIF

! effective cloudiness accounting for condensed water
        ZCLDLD(JL,JK) = PCLFR(JL,IKL)*(1.0_JPRB-EXP(-ZMSALD*ZZFLWP-ZMSAID*&
         & ZZFIWP))  
        ZCLDLU(JL,JK) = PCLFR(JL,IKL)*(1.0_JPRB-EXP(-ZMSALU*ZZFLWP-ZMSAIU*&
         & ZZFIWP))  
      ENDIF
    ENDDO

  ELSE

!          2.5    CLOUD LONGWAVE OPTICAL PROPERTIES FOR RRTM
!                 ------------------------------------------

!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Savijarvi (1998)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

! No need for a fixed diffusivity factor, accounted for spectrally below
! The detailed spectral structure does not require defining upward and
! downward effective optical properties
    ICE_WATER_CLOUD_TOT_EMIS = -1
    IF (NLWLIQOPT == 0 .OR. NLWLIQOPT >= 3 ) THEN
      ! water cloud total emissivity a la Smith and Shi (1992)
      IF (NLWICEOPT == 0) THEN
        ! ice cloud spectral emissivity a la Smith & Shi (1992)
        ICE_WATER_CLOUD_TOT_EMIS = 1
      ELSEIF (NLWICEOPT == 1) THEN
        ! ice cloud spectral emissivity a la Ebert-Curry (1992)
        ICE_WATER_CLOUD_TOT_EMIS = 2
      ELSEIF (NLWICEOPT == 2) THEN
        ! ice cloud spectral emissivity a la Fu & Liou (1993)
        ICE_WATER_CLOUD_TOT_EMIS = 3
      ELSEIF (NLWICEOPT == 3) THEN
        ! ice cloud spectral emissivity a la Fu et al (1998) including 
        ! parametrisation for LW scattering effect  
        ICE_WATER_CLOUD_TOT_EMIS = 4
      ENDIF    
    ELSEIF (NLWLIQOPT == 1) THEN
      ! water cloud spectral emissivity a la Savijarvi (1997)
      IF (NLWICEOPT == 0) THEN
        ! ice cloud spectral emissivity a la Smith & Shi (1992)
        ICE_WATER_CLOUD_TOT_EMIS = 5
      ELSEIF (NLWICEOPT == 1) THEN
        ! ice cloud spectral emissivity a la Ebert-Curry (1992)
        ICE_WATER_CLOUD_TOT_EMIS = 6
      ELSEIF (NLWICEOPT == 2) THEN
        ! ice cloud spectral emissivity a la Fu & Liou (1993)
        ICE_WATER_CLOUD_TOT_EMIS = 7
      ELSEIF (NLWICEOPT == 3) THEN
        ! ice cloud spectral emissivity a la Fu et al (1998) including 
        ! parametrisation for LW scattering effect  
        ICE_WATER_CLOUD_TOT_EMIS = 8
      ENDIF    
    ELSEIF (NLWLIQOPT == 2) THEN
      ! water cloud spectral emissivity a la Lindner and Li (2000)
      IF (NLWICEOPT == 0) THEN
        ! ice cloud spectral emissivity a la Smith & Shi (1992)
        ICE_WATER_CLOUD_TOT_EMIS = 9
      ELSEIF (NLWICEOPT == 1) THEN
        ! ice cloud spectral emissivity a la Ebert-Curry (1992)
        ICE_WATER_CLOUD_TOT_EMIS = 10
      ELSEIF (NLWICEOPT == 2) THEN
        ! ice cloud spectral emissivity a la Fu & Liou (1993)
        ICE_WATER_CLOUD_TOT_EMIS = 11
      ELSEIF (NLWICEOPT == 3) THEN
        ! ice cloud spectral emissivity a la Fu et al (1998) including 
        ! parametrisation for LW scattering effect  
        ICE_WATER_CLOUD_TOT_EMIS = 12
      ENDIF    
    ENDIF  
    DO JRTM=1,16
      DO JL = KIDIA,KFDIA
        ZTAUCLD(JL,JK,JRTM) = 0.0_JPRB
      ENDDO
    ENDDO
    SELECT CASE(ICE_WATER_CLOUD_TOT_EMIS)
    CASE (1)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
      
            ! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
            ZRSALD= 0.144_JPRB*ZMULTL / 1.66_JPRB
            
            ! ice cloud spectral emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
            ZRSAID= 0.103_JPRB*ZMULTI / 1.66_JPRB
            
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (2)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
            ZRSALD= 0.144_JPRB*ZMULTL / 1.66_JPRB
            
            ! ice cloud spectral emissivity a la Ebert-Curry (1992)
            ZRSAID= REBCUH(JRTM)+REBCUG(JRTM)/ZRADIP(JL)
            
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (3)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
            ZRSALD= 0.144_JPRB*ZMULTL / 1.66_JPRB
            
            ! ice cloud spectral emissivity a la Fu & Liou (1993)
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZRSAID = RFULIO(JRTM,1) + Z1RADI&
            & *(RFULIO(JRTM,2) + Z1RADI * RFULIO(JRTM,3))  
             
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (4)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
            ZRSALD= 0.144_JPRB*ZMULTL / 1.66_JPRB
            
            ! ice cloud spectral emissivity a la Fu et al (1998) including 
            ! parametrisation for LW scattering effect  
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZRSAIE = RFUETA(JRTM,1) + Z1RADI&
             &*(RFUETA(JRTM,2) + Z1RADI * RFUETA(JRTM,3)) 
            ZRSAIA = Z1RADI*(RFUETB(JRTM,1)&
             & +ZDESR(JL)*( RFUETB(JRTM,2)&
             & +ZDESR(JL)*( RFUETB(JRTM,3) +ZDESR(JL)* RFUETB(JRTM,4))))
            ZRSAIG = RFUETC(JRTM,1) +ZDESR(JL)*( RFUETC(JRTM,2)&
             & +ZDESR(JL)*( RFUETC(JRTM,3) +ZDESR(JL)* RFUETC(JRTM,4))) 
            ZRSAIF = 0.5_JPRB + ZRSAIG*( 0.3738_JPRB&
             & + ZRSAIG*( 0.0076_JPRB + ZRSAIG*0.1186_JPRB ) )
!            ZRSAID = (1.0_JPRB - ZRSAIA/ZRSAIE * ZRSAIF) * ZRSAIE
            ZRSAID = ZRSAIE - ZRSAIA*ZRSAIF
         
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (5)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Savijarvi (1997)
            ZRSALD= RHSAVI(JRTM,1) + ZRADLP(JL)&
            & *(RHSAVI(JRTM,2) + ZRADLP(JL)*RHSAVI(JRTM,3))  
             
            ! ice cloud spectral emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
            ZRSAID= 0.103_JPRB*ZMULTI / 1.66_JPRB
            
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (6)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Savijarvi (1997)
            ZRSALD= RHSAVI(JRTM,1) + ZRADLP(JL)&
            & *(RHSAVI(JRTM,2) + ZRADLP(JL)*RHSAVI(JRTM,3))  
             
            ! ice cloud spectral emissivity a la Ebert-Curry (1992)
            ZRSAID= REBCUH(JRTM)+REBCUG(JRTM)/ZRADIP(JL)
            
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (7)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Savijarvi (1997)
            ZRSALD= RHSAVI(JRTM,1) + ZRADLP(JL)&
            & *(RHSAVI(JRTM,2) + ZRADLP(JL)*RHSAVI(JRTM,3))  
             
            ! ice cloud spectral emissivity a la Fu & Liou (1993)
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZRSAID = RFULIO(JRTM,1) + Z1RADI&
            & *(RFULIO(JRTM,2) + Z1RADI * RFULIO(JRTM,3))  
             
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (8)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Savijarvi (1997)
            ZRSALD= RHSAVI(JRTM,1) + ZRADLP(JL)&
            & *(RHSAVI(JRTM,2) + ZRADLP(JL)*RHSAVI(JRTM,3))  
             
            ! ice cloud spectral emissivity a la Fu et al (1998) including 
            ! parametrisation for LW scattering effect  
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZRSAIE = RFUETA(JRTM,1) + Z1RADI&
            &*(RFUETA(JRTM,2) + Z1RADI * RFUETA(JRTM,3)) 
            ZRSAIA = Z1RADI*(RFUETB(JRTM,1) +ZDESR(JL)*( RFUETB(JRTM,2) +&
                                            &ZDESR(JL)*( RFUETB(JRTM,3) +&
                                            &ZDESR(JL)* RFUETB(JRTM,4))))
            ZRSAIG = RFUETC(JRTM,1) +ZDESR(JL)*( RFUETC(JRTM,2) +&
                                    &ZDESR(JL)*( RFUETC(JRTM,3) +&
                                    &ZDESR(JL)* RFUETC(JRTM,4))) 
            ZRSAIF = 0.5_JPRB + ZRSAIG*( 0.3738_JPRB + ZRSAIG*( 0.0076_JPRB + ZRSAIG*0.1186_JPRB ) )
            ZRSAID = (1.0_JPRB - ZRSAIA/ZRSAIE * ZRSAIF) * ZRSAIE
         
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (9)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = 1.0_JPRB / ZRADLP(JL)
            ZEXTCF = RLILIA(JRTM,1)+ZRADLP(JL)*RLILIA(JRTM,2)+ Z1RADL*&
            & (RLILIA(JRTM,3) + Z1RADL*(RLILIA(JRTM,4) + Z1RADL*&
            & RLILIA(JRTM,5) ))  
            Z1MOMG = RLILIB(JRTM,1) + Z1RADL*RLILIB(JRTM,2)&
            & + ZRADLP(JL) *(RLILIB(JRTM,3) + ZRADLP(JL)*RLILIB(JRTM,4) )
            ZRSALD = Z1MOMG * ZEXTCF
         
            ! ice cloud spectral emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
            ZRSAID= 0.103_JPRB*ZMULTI / 1.66_JPRB
            
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (10)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = 1.0_JPRB / ZRADLP(JL)
            ZEXTCF = RLILIA(JRTM,1)+ZRADLP(JL)*RLILIA(JRTM,2)+ Z1RADL*&
            & (RLILIA(JRTM,3) + Z1RADL*(RLILIA(JRTM,4) + Z1RADL*&
            & RLILIA(JRTM,5) ))  
            Z1MOMG = RLILIB(JRTM,1) + Z1RADL*RLILIB(JRTM,2)&
            & + ZRADLP(JL) *(RLILIB(JRTM,3) + ZRADLP(JL)*RLILIB(JRTM,4) )
            ZRSALD = Z1MOMG * ZEXTCF
         
            ! ice cloud spectral emissivity a la Ebert-Curry (1992)
            ZRSAID= REBCUH(JRTM)+REBCUG(JRTM)/ZRADIP(JL)

            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)            
          ENDIF
        ENDDO
      ENDDO
    CASE (11)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = 1.0_JPRB / ZRADLP(JL)
            ZEXTCF = RLILIA(JRTM,1)+ZRADLP(JL)*RLILIA(JRTM,2)+ Z1RADL*&
            & (RLILIA(JRTM,3) + Z1RADL*(RLILIA(JRTM,4) + Z1RADL*&
            & RLILIA(JRTM,5) ))  
            Z1MOMG = RLILIB(JRTM,1) + Z1RADL*RLILIB(JRTM,2)&
            & + ZRADLP(JL) *(RLILIB(JRTM,3) + ZRADLP(JL)*RLILIB(JRTM,4) )
            ZRSALD = Z1MOMG * ZEXTCF
         
            ! ice cloud spectral emissivity a la Fu & Liou (1993)
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZRSAID = RFULIO(JRTM,1) + Z1RADI&
            & *(RFULIO(JRTM,2) + Z1RADI * RFULIO(JRTM,3))  
             
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE (12)
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
    
            ! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = 1.0_JPRB / ZRADLP(JL)
            ZEXTCF = RLILIA(JRTM,1)+ZRADLP(JL)*RLILIA(JRTM,2)+ Z1RADL*&
            & (RLILIA(JRTM,3) + Z1RADL*(RLILIA(JRTM,4) + Z1RADL*&
            & RLILIA(JRTM,5) ))  
            Z1MOMG = RLILIB(JRTM,1) + Z1RADL*RLILIB(JRTM,2)&
            & + ZRADLP(JL) *(RLILIB(JRTM,3) + ZRADLP(JL)*RLILIB(JRTM,4) )
            ZRSALD = Z1MOMG * ZEXTCF
         
            ! ice cloud spectral emissivity a la Fu et al (1998) including 
            ! parametrisation for LW scattering effect  
            Z1RADI = 1.0_JPRB / ZDESR(JL)
            ZRSAIE = RFUETA(JRTM,1) + Z1RADI&
             &*(RFUETA(JRTM,2) + Z1RADI * RFUETA(JRTM,3)) 
            ZRSAIA = Z1RADI*(RFUETB(JRTM,1) +ZDESR(JL)*( RFUETB(JRTM,2) +&
                                            &ZDESR(JL)*( RFUETB(JRTM,3) +&
                                            &ZDESR(JL)* RFUETB(JRTM,4))))
            ZRSAIG = RFUETC(JRTM,1) +ZDESR(JL)*( RFUETC(JRTM,2) +&
                                    &ZDESR(JL)*( RFUETC(JRTM,3) +&
                                    &ZDESR(JL)* RFUETC(JRTM,4))) 
            ZRSAIF = 0.5_JPRB + ZRSAIG*( 0.3738_JPRB + ZRSAIG*( 0.0076_JPRB + ZRSAIG*0.1186_JPRB ) )
            ZRSAID = (1.0_JPRB - ZRSAIA/ZRSAIE * ZRSAIF) * ZRSAIE
         
            ZTAUCLD(JL,JK,JRTM) = ZRSALD*ZFLWP(JL)+ZRSAID*ZFIWP(JL)

          ENDIF
        ENDDO
      ENDDO
    CASE DEFAULT
      WRITE(NULOUT,*) ' NLWLIQOPT =',NLWLIQOPT,' NLWICEOPT =',NLWICEOPT
      CALL ABOR1("radlsw : incorrect value for emissivity")
    END SELECT
        
    ! Diffusivity correction within clouds a la Savijarvi
    IF (LDIFFC) THEN
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
            ZDIFFD=MIN(MAX(1.517_JPRB-0.156_JPRB*LOG(ZTAUCLD(JL,JK,JRTM)) , 1.0_JPRB), 2.0_JPRB)
            ZTAUCLD(JL,JK,JRTM) = ZTAUCLD(JL,JK,JRTM)*ZDIFFD
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO JRTM=1,16
        DO JL = KIDIA,KFDIA
          IF (ZFLWP(JL)+ZFIWP(JL) > REPSCW) THEN
            ZDIFFD=1.66_JPRB
            ZTAUCLD(JL,JK,JRTM) = ZTAUCLD(JL,JK,JRTM)*ZDIFFD
          ENDIF
        ENDDO
      ENDDO
    ENDIF
        

  ENDIF

ENDDO

NUAER = NUA
NTRAER = NTRA

!     ------------------------------------------------------------------

!          2.6    SCALING OF OPTICAL THICKNESS
!                 SPECTRALLY, ACCOUNTING FOR VERTICAL VARIABILITY

JEXPLR=NLAYINH
JXPLDN=2*JEXPLR+1

IF (NINHOM == 1) THEN
!-- simple scaling a la Tiedtke (1996) with RSWINHF in SW and RLWINHF in LW
  DO JSW=1,NSW
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZTAU(JL,JSW,JK)=ZTAU(JL,JSW,JK) * RSWINHF
      ENDDO
    ENDDO
  ENDDO

  DO JRTM=1,16
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZTAUCLD(JL,JK,JRTM)=ZTAUCLD(JL,JK,JRTM) * RLWINHF
      ENDDO
    ENDDO
  ENDDO

ELSEIF (JEXPLR /= 0) THEN
  DO JSW=1,NSW
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZSQUAR(JL,JK)=0.0_JPRB
        ZVARIA(JL,JK)=1.0_JPRB
      ENDDO
    ENDDO
!-- range should be defined from Hogan & Illingworth
    DO JK=1+JEXPLR,KLEV-JEXPLR
      DO JL=KIDIA,KFDIA
!        ZAVDP(JL)=0.0_JPRB
        ZAVTO(JL)=0.0_JPRB
        ZSQTO(JL)=0.0_JPRB
      ENDDO
      DO JKI=JK-JEXPLR,JK+JEXPLR
        IKI=KLEV+1-JKI
        DO JL=KIDIA,KFDIA
!          ZAVDP(JL)=ZAVDP(JL)+PDP(JL,IKI)/RG
          ZAVTO(JL)=ZAVTO(JL)+ZTAU(JL,JSW,JKI)
        ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
!        ZAVTO(JL)=ZAVTO(JL)/ZAVDP(JL)
        ZAVTO(JL)=ZAVTO(JL)/JXPLDN
      ENDDO
      DO JKI=JK-JEXPLR,JK+JEXPLR
        IKI=KLEV+1-JKI
        DO JL=KIDIA,KFDIA
!          ZSQTO(JL)=ZSQTO(JL)+(ZTAU(JL,JSW,JKI)/PDP(JL,IKI)-ZAVTO(JL))**2
          ZSQTO(JL)=ZSQTO(JL)+(ZTAU(JL,JSW,JKI)-ZAVTO(JL))**2
        ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
        ZSQTO(JL)=SQRT(ZSQTO(JL)/(JXPLDN*(JXPLDN-1)))
        IF (ZAVTO(JL) > 0.0_JPRB) THEN
          ZVARIA(JL,JK)=(ZSQTO(JL)/ZAVTO(JL))**2
          ZSQUAR(JL,JK)=EXP(-ZVARIA(JL,JK))
        ELSE
          ZVARIA(JL,JK)=0.0_JPRB
          ZSQUAR(JL,JK)=1.0_JPRB
        ENDIF

!-- scaling a la Barker
        IF (NINHOM ==2) THEN
          ZTAU(JL,JSW,JK)=ZTAU(JL,JSW,JK)*ZSQUAR(JL,JK)

!-- scaling a la Cairns et al.
        ELSEIF (NINHOM == 3) THEN
          ZVI=ZVARIA(JL,JK) 
          ZTAU(JL,JSW,JK)  = ZTAU(JL,JSW,JK)/(1.0_JPRB+ZVI)
          ZOMEGA(JL,JSW,JK)= ZOMEGA(JL,JSW,JK)&
            &   /(1.0_JPRB + ZVI*(1.0_JPRB-ZOMEGA(JL,JSW,JK) ) )
          ZCG(JL,JSW,JK)   = ZCG(JL,JSW,JK)&
            & *(1.0_JPRB+ZVI*(1.0_JPRB-ZOMEGA(JL,JSW,JK)))&
            & /(1.0_JPRB+ZVI*(1.0_JPRB-ZOMEGA(JL,JSW,JK)*ZCG(JL,JSW,JK)))
        ENDIF
      ENDDO
!      JL=KIDIA
!      print 9261,JSW,JK,ZTAU(JL,JSW,JK),ZAVTO(JL),ZSQTO(JL),ZVARIA(JL,JK),ZSQUAR(JL,JK)
9261   FORMAT(1X,'Varia1 ',2I3,7F10.4)
    ENDDO
  ENDDO


  DO JRTM=1,16
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZSQUAR(JL,JK)=0.0_JPRB
        ZVARIA(JL,JK)=1.0_JPRB
      ENDDO
    ENDDO
!-- range to be defined from Hogan & Illingworth
    DO JK=1+JEXPLR,KLEV-JEXPLR
      DO JL=KIDIA,KFDIA
!        ZAVDP(JL)=0.0_JPRB
        ZAVTO(JL)=0.0_JPRB
        ZSQTO(JL)=0.0_JPRB
      ENDDO
      DO JKI=JK-JEXPLR,JK+JEXPLR
        IKI=KLEV+1-JKI
        DO JL=KIDIA,KFDIA
!          ZAVDP(JL)=ZAVDP(JL)+PDP(JL,IKI)/RG
          ZAVTO(JL)=ZAVTO(JL)+ZTAUCLD(JL,JKI,JRTM)
        ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
!        ZAVTO(JL)=ZAVTO(JL)/ZAVDP(JL)
        ZAVTO(JL)=ZAVTO(JL)/JXPLDN
      ENDDO
      DO JKI=JK-JEXPLR,JK+JEXPLR
        IKI=KLEV+1-JKI
        DO JL=KIDIA,KFDIA
!          ZSQTO(JL)=ZSQTO(JL)+(ZTAUCLD(JL,JKI,JRTM)/PDP(JL,IKI)-ZAVTO(JL))**2
            ZSQTO(JL)=ZSQTO(JL)+(ZTAUCLD(JL,JKI,JRTM)-ZAVTO(JL))**2
        ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
        ZSQTO(JL)=SQRT(ZSQTO(JL)/(JXPLDN*(JXPLDN-1)))
        IF (ZAVTO(JL) > 0.0_JPRB) THEN
          ZVARIA(JL,JK)=(ZSQTO(JL)/ZAVTO(JL))**2
          ZSQUAR(JL,JK)=EXP(-ZVARIA(JL,JK))
        ELSE
          ZVARIA(JL,JK)=0.0_JPRB
          ZSQUAR(JL,JK)=1.0_JPRB
        ENDIF

!-- scaling a la Barker
        IF (NINHOM ==2) THEN
          ZTAUCLD(JL,JK,JRTM)=ZTAUCLD(JL,JK,JRTM)*ZSQUAR(JL,JK)

!-- scaling a la Cairns et al.
        ELSEIF (NINHOM == 3) THEN
          ZVI=ZVARIA(JL,JK) 
          ZTAUCLD(JL,JK,JRTM)=ZTAUCLD(JL,JK,JRTM)/(1.0_JPRB+ZVI)
        ENDIF
      ENDDO
!      JL=KIDIA
!      print 9262,JRTM,JK,ZTAUCLD(JL,JK,JRTM),ZAVTO(JL),ZSQTO(JL),ZVARIA(JL,JK),ZSQUAR(JL,JK)
9262   FORMAT(1X,'Varia2 ',2I3,7F10.4)
    ENDDO
  ENDDO
ENDIF



!     ------------------------------------------------------------------

!*         2.7    DIFFUSIVITY FACTOR OR SATELLITE VIEWING ANGLE
!                 ---------------------------------------------

DO JL = KIDIA,KFDIA
  ZVIEW(JL) = DIFF
ENDDO

!     ------------------------------------------------------------------

!*         3.     CALL LONGWAVE RADIATION CODE
!                 ----------------------------

!*         3.1    FULL LONGWAVE RADIATION COMPUTATIONS
!                 ------------------------------------

IF ( .NOT. LRRTM .OR. LPHYLIN) THEN
!* Computation of the cloud top level up to which cloud contribution
!  is computed in lwc.F90
  IF (LPHYLIN .AND. LWLCLHR) THEN

    DO JL = KIDIA,KFDIA
      ILEV(JL)=KLEV
    ENDDO
    DO JK = 1 , KLEV
      DO JL = KIDIA,KFDIA
        ZCOND(JL,JK) = PQLWP(JL,JK)+PQIWP(JL,JK)
        IF (PCLFR(JL,JK)>=1.E-04_JPRB&
         & .AND. ZCOND(JL,JK)>=0.5E-08_JPRB) THEN
          ILEV(JL) = MIN(ILEV(JL),JK)
        ENDIF
      ENDDO
    ENDDO

    ILEVR=KLEV
    DO JL = KIDIA,KFDIA
      ILEVR = MAX(NLEVLWC,MIN(ILEVR,ILEV(JL)))
    ENDDO
!orig    ITOPC = KLEV-ILEVR
    ITOPC = MIN(KLEV,KLEV-ILEVR+7)
  ELSE
    ITOPC = KLEV
  ENDIF
!* End of the cloud top computation

!or   IF (LPHYLIN .AND. LH2OCO2) THEN
  IF (LH2OCO2) THEN
    IUAER  = 9
  ELSE
    IUAER  = NUAER
  ENDIF

  CALL LW&
   & ( YDML_PHY_RAD,YDEPHLI,YDPHNC, KIDIA , KFDIA , KLON  , KLEV , KMODE, ITOPC,IUAER,&
   & PCCO2 , PCLFR , ZCLDLD, ZCLDLU,&
   & PDP   , ZDT0  , ZEMIS , ZEMIW,&
   & ZPMB  , POZON , ZTL,&
   & PAER  , ZTAVE , ZVIEW , PQ,&
   & ZEMIT , ZFLUX , ZFLUC&
   & )  

ELSE

!*         3.2    FULL LONGWAVE RADIATION COMPUTATIONS - RRTM
!                 ------------------------------------   ----

!  i)  pass ZOZN (ozone mass mixing ratio) to RRTM; remove pressure
!      weighting applied to POZON in driverMC (below)
!  ii) pass ZEMIS and ZEMIW to RRTM; return ZEMIT from RRTM
!  iii)pass ZTAUCLD, cloud optical depths (water+ice) to RRTM, 
!      computed from equations above
!  iv) pass ECRT arrays to RRTM arrays in interface routine ECRTATM
!      in module rrtm_ecrt.f

  DO JL = KIDIA,KFDIA
    DO JK = 1, KLEV
      ZOZN(JL,JK) = POZON(JL,JK)/PDP(JL,JK)
      ZCO2(JL,JK) = PCCO2
      ZCH4(JL,JK) = RCH4
      ZN2O(JL,JK) = RN2O
      ZNO2(JL,JK) = 0._JPRB
      ZC11(JL,JK) = RCFC11
      ZC12(JL,JK) = RCFC12
      ZC22(JL,JK) = RCFC22
      ZCL4(JL,JK) = RCCL4
    ENDDO
  ENDDO

  CALL RRTM_RRTM_140GP&
   & (YDDIMV, YDML_PHY_RAD%YRERAD, KIDIA , KFDIA , KLON  , KLEV,&
   & YDML_PHY_RAD%YRESWRT%RAER, PAER  , PAPH , PAP  ,&
   & PTS   , PTH  , PT   ,&
   & ZEMIS , ZEMIW,&
   & PQ    , ZCO2 , ZCH4 , ZN2O   , ZNO2, ZC11, ZC12, ZC22, ZCL4, ZOZN, ZCLDSW, ZTAUCLD,&
   & ZEMIT , ZFLUX, ZFLUC, ZTCLEAR&
   & )  
ENDIF

!     ------------------------------------------------------------------

!*         4.     CALL SHORTWAVE RADIATION CODE
!                 -----------------------------

ZRMUZ=0.0_JPRB
DO JL = KIDIA,KFDIA
  ZRMUZ = MAX (ZRMUZ, ZMU0(JL))
ENDDO

IF (NSTEP == 0 .AND. LEDBUG .AND. ZMU0(KIDIA) > 0.0_JPRB) THEN
  WRITE(NULOUT,'(4E15.8)') PRII0,PCCO2,ZPSOL(KIDIA),ZMU0(KIDIA)
  WRITE(NULOUT,'("ZALBD ",6E15.8)') (ZALBD(KIDIA,JSW),JSW=1,NSW)
  WRITE(NULOUT,'("ZALBP ",6E15.8)') (ZALBP(KIDIA,JSW),JSW=1,NSW)
  WRITE(NULOUT,'("PQ    ",10E12.5)') (PQ(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("PQS   ",10E12.5)') (PQS(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("PDP   ",10E12.5)') (PDP(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZPMB  ",10E12.5)') (ZPMB(KIDIA,JK),JK=1,KLEV+1)
  WRITE(NULOUT,'("ZTAVE ",10E12.5)') (ZTAVE(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZCLDSW",10E12.5)') (ZCLDSW(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZTAU  ",10E12.5)') ((ZTAU(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
  WRITE(NULOUT,'("ZCG   ",10E12.5)') ((ZCG(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
  WRITE(NULOUT,'("ZOMEGA",10E12.5)') ((ZOMEGA(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
  WRITE(NULOUT,'("ZOZ   ",10E12.5)') (ZOZ(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("PAER  ",10E12.5)') ((PAER(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
ENDIF

IF (NSTEP == 0 .AND. LEDBUG .AND. ZMU0(KIDIA) > 0.0_JPRB) THEN
  WRITE(NULOUT,'(4E15.8)') PRII0,PCCO2,ZPSOL(KIDIA),ZMU0(KIDIA)
  WRITE(NULOUT,'("ZALBD ",6E15.8)') (ZALBD(KIDIA,JSW),JSW=1,NSW)
  WRITE(NULOUT,'("ZALBP ",6E15.8)') (ZALBP(KIDIA,JSW),JSW=1,NSW)
  WRITE(NULOUT,'("PQ    ",10E12.5)') (PQ(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("PQS   ",10E12.5)') (PQS(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("PDP   ",10E12.5)') (PDP(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZPMB  ",10E12.5)') (ZPMB(KIDIA,JK),JK=1,KLEV+1)
  WRITE(NULOUT,'("ZTAVE ",10E12.5)') (ZTAVE(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZCLDSW",10E12.5)') (ZCLDSW(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZTAU  ",10E12.5)') ((ZTAU(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
  WRITE(NULOUT,'("ZCG   ",10E12.5)') ((ZCG(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
  WRITE(NULOUT,'("ZOMEGA",10E12.5)') ((ZOMEGA(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
  WRITE(NULOUT,'("ZOZ   ",10E12.5)') (ZOZ(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("PAER  ",10E12.5)') ((PAER(KIDIA,JSW,JK),JK=1,KLEV),JSW=1,NSW)
ENDIF
CALL SW&
 & ( YDML_PHY_RAD,YDPHNC,YDECLD, KIDIA , KFDIA , KLON  , KLEV  , KAER,&
 & PRII0 , PCCO2 , ZPSOL , ZALBD , ZALBP , PQ   , PQS,&
 & ZMU0  , ZCG   , ZCLDSW, ZOMEGA, ZOZ  , ZPMB,&
 & ZTAU  , ZTAVE , PAER,&
 & ZFSDWN, ZFSUP , ZFCDWN, ZFCUP,&
 & ZFSDNN, ZFSDNV, ZFSUPN, ZFSUPV,&
 & ZFCDNN, ZFCDNV, ZFCUPN, ZFCUPV,&
 & ZSUDU , ZUVDF , ZPARF ,ZPARCF, ZDIFFS, ZDIRFS,&
 & LDDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST&
   & )
PFSDNV=ZFSDNV
PFSDNN=ZFSDNN
IF (SIZE(PSFSWDIR,2)>1) THEN
  PSFSWDIR= ZDIRFS
  PSFSWDIF= ZDIFFS
ELSE
  PSFSWDIR (:,1) = ZFSDNV(:) + ZFSDNN(:)
  PSFSWDIF (:,:) = 0.
ENDIF

IF (NSTEP == 0 .AND. LEDBUG .AND. ZMU0(KIDIA) > 0.0_JPRB) THEN
  WRITE(NULOUT,'("ZFSDWN",10E12.5)') (ZFSDWN(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZFSUP ",10E12.5)') (ZFSUP (KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZFCDWN",10E12.5)') (ZFCDWN(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZFCUP ",10E12.5)') (ZFCUP (KIDIA,JK),JK=1,KLEV)
  LEDBUG=.FALSE. 
ENDIF
IF (NSTEP == 0 .AND. LEDBUG .AND. ZMU0(KIDIA) > 0.0_JPRB) THEN
  WRITE(NULOUT,'("ZFSDWN",10E12.5)') (ZFSDWN(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZFSUP ",10E12.5)') (ZFSUP (KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZFCDWN",10E12.5)') (ZFCDWN(KIDIA,JK),JK=1,KLEV)
  WRITE(NULOUT,'("ZFCUP ",10E12.5)') (ZFCUP (KIDIA,JK),JK=1,KLEV)
  LEDBUG=.FALSE. 
ENDIF
!     ------------------------------------------------------------------

!*         5.     FILL UP THE MODEL NET LW AND SW RADIATIVE FLUXES
!                 ------------------------------------------------

DO JKL = 1 , KLEV+1
  JK = KLEV+1 + 1 - JKL
  DO JL = KIDIA,KFDIA
    PFLS(JL,JKL) = ZFSDWN(JL,JK) - ZFSUP(JL,JK)
    PFLT(JL,JKL) = - ZFLUX(JL,1,JK) - ZFLUX(JL,2,JK)
    PFCS(JL,JKL) = ZFCDWN(JL,JK) - ZFCUP(JL,JK)
    PFCT(JL,JKL) = - ZFLUC(JL,1,JK) - ZFLUC(JL,2,JK)
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  PFRSOD(JL)=ZFSDWN(JL,1)
  PEMIT (JL)=ZEMIT (JL)
  PSUDU (JL)=ZSUDU (JL)
  PUVDF (JL)=ZUVDF (JL)
  PPARF (JL)=ZPARF (JL)
  PPARCF(JL)=ZPARCF(JL)
  PTINCF(JL)=PRII0 * ZMU0(JL) 
ENDDO
!print 9501,(PUVDF(JL),JL=KIDIA,KFDIA)
9501 FORMAT(1X,'Radlsw PUVDF: ',30F6.1)
!print 9502,(PPARF(JL),JL=KIDIA,KFDIA)
9502 FORMAT(1X,'Radlsw PPARF: ',30F6.1)

!     --------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADLSW',1,ZHOOK_HANDLE)
END SUBROUTINE RADLSW
