! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SRC &
 &( YDGEOMETRY, YDMODEL,KIDIA, KFDIA, KLON , KTDIA, KLEV, KTILES, KSTART, KSTEP ,KSTGLO,&
 &  KTRAC, KAERO, &
 &  PAPHI , &
 &  PAERGUST, PALTH ,&
 &  PDMS,&
 &  PAPH , PAP  , PCI  , PGELAM, PGELAT, PGEMU , PFRTI , &
 &  PLSM , PCLK , PQ   , PRHO , PSNS  , PT    , PTS   , PTSPHY, KCHEM,&
 &  PWIND, PWS1 , PSST, &
 &  PDMSO, PLDAY, PLISS, PSO2  , PTDMS, PAERFLX, PCFLX , PCEN  , PTENC,&
 &  PODMS)

!*** * AER_SRC* - SOURCE TERMS FOR TEST AEROSOLS

!**   INTERFACE.
!     ----------
!          *AER_SRC* IS CALLED FROM *AER_PHY2*


!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        ORIGINAL : 2004-04-27

!     MODIFICATIONS.
!     --------------
!        JJMorcrette 20110707  Revision for handling volcanic aerosol emission
!        JJMorcrette 20111222  Revision of dust emission formulation
!        JJMorcrette 20130619  oceanic DMS
!        K. Yessad (July 2014): Move some variables.
!        S.Rémy      20150126  Biomass burning emissions injection heights
!        S.Rémy      20150126  MACCity SO2 emissions : no need to renormalize with mass
!        S.Rémy      20150130  Add P.Nabat's scheme for dust emissions
!                              (Marticorena et al., 1995, Kok et al., 2011)
!        S.Rémy      20150609  SOA Emissions scaled on CO emissions (Spracklen
!        et al., 2011)
!        S.Rémy      20160204  Code cleaning
!        S.Rémy      20170420  Add Grythe SS emissions from P.Nabat
!        S.Rémy      20170420  SOA emissions = CO scaled + EDGAR
!        S.Rémy      20170421  Lifted SO2 emissions above volcanoes
!        S.Rémy      20170512  Diurnal cycle for BB aerosol emissions +
!        externalise diurnal cycle
!        S.Rémy      201707    New dust scheme, LAERDUSTSOURCE and
!        LAERDUST_NEWBIN
!        R.Hogan     20190114  Change LE4ALB to NALBEDOSCHEME
!        R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!-----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMCST       , ONLY : YRCST
USE YOETHF       , ONLY : YRTHF
USE YOMRIP0      , ONLY : NINDAT, NSSSSS
IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON, KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV, KTDIA, KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP, KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV), PALTH(KLON,0:KLEV) 
REAL(KIND=JPRB),   INTENT(IN)    :: PAERGUST(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV), PAPH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON), PGELAT(KLON), PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLON), PLSM(KLON), PCLK(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV), PRHO(KLON,KLEV), PSNS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS1(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAERFLX(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFLX(KLON,KTRAC)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMSO(KLON), PLDAY(KLON), PLISS(KLON), PSO2(KLON), PTDMS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PODMS(KLON)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JAER, JK, JL, JN, JNT, JNM
INTEGER(KIND=JPIM) :: IFLAG, INBAER, IHTST, ITEST, IWAVL
INTEGER(KIND=JPIM) :: JS

LOGICAL ::  LLDUST(KLON), LLPDUSTS(KLON), LLDDACT(KLON)
!LOGICAL :: LLDVOL(KLON)

REAL(KIND=JPRB) :: ZFAERO(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)     , ZAEROCLIS(KLON,KLEV,2) 
REAL(KIND=JPRB) :: ZAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO), ZTAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)

REAL(KIND=JPRB) :: ZGLAT(KLON), ZGLON(KLON)
REAL(KIND=JPRB) :: ZHDD, ZHSS
REAL(KIND=JPRB) :: ZDETAH(KLON,KLEV), ZETA(KLON,KLEV) , ZETAH(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV) , ZRHCL(KLON,KLEV), ZTH(KLON,0:KLEV)

REAL(KIND=JPRB) :: ZFLX_SDUST(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
REAL(KIND=JPRB) :: ZFLX_SSALT(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(1))
REAL(KIND=JPRB) :: ZSCC2(KLON) , ZDEP2(KLON) , ZLTS2(KLON), ZLTSMIN(KLON), ZLTSMAX(KLON)
REAL(KIND=JPRB) :: ZWNDDU(KLON), ZWNDSS(KLON), ZWND3(KLON) 
REAL(KIND=JPRB) :: ZDUEMPOT(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
REAL(KIND=JPRB) :: ZDEGRAD , ZFSWET , ZSWETN, ZSWETN2
REAL(KIND=JPRB) :: ZSO2MSS
REAL(KIND=JPRB) :: ZEPSISS, ZEPSIDD, ZEPSIRA, ZEPSSNO, ZEPSARE
!REAL(KIND=JPRB) :: ZHVO

!-- volcano-related variables
INTEGER(KIND=JPIM) :: INDLAT(KLON)
INTEGER(KIND=JPIM) :: IBASPL, ITOPPL, IGLGLO, IGVOLC, IGVOLE
INTEGER(KIND=JPIM) :: JVOLC, JVOLE
INTEGER(KIND=JPIM) :: IVDATES, IVDAYS, IVSECNDS, IVDAYE, IVSECNDE
INTEGER(KIND=JPIM) :: IYY, IMM, IDD, IMDATE 
INTEGER(KIND=JPIM) :: IY0, IM0, ID0, INC, IMON(12)

REAL(KIND=JPRB)    :: ZDIST(KLON), ZGRDLON(KLON) , ZGELAT(KLON), ZGDLAT(KLON), ZGDLON(KLON)
REAL(KIND=JPRB)    :: ZVOLCASH(KLON), ZVOLEASH(KLON)
REAL(KIND=JPRB)    :: ZVOLCSO2(KLON), ZVOLESO2(KLON)
REAL(KIND=JPRB)    :: ZAREA, ZCIRCLE, ZDEEPLUME, ZDLAT, ZDLON, ZEQUATOR, ZLATSQ, ZLONSQ
REAL(KIND=JPRB)    :: ZGELAA, ZGELAB, ZGELAV, ZGELOA, ZGELOB, ZGELOV, ZPIH, ZPI2, ZVOLUME
REAL(KIND=JPRB)    :: ZGRDLAT, ZGRDLAT2, ZGRDLON2, ZINCLAT, ZA, ZC, Z1GP

REAL(KIND=JPRB)    :: ZLONE, ZLONGB, ZLONW, ZLAT
REAL(KIND=JPRB)    :: ZREFRAD, ZREFSPD, ZRADREF

!-- QnD oceanic DMS
REAL(KIND=JPRB)    :: ZCOS0, ZSIN0, ZRAD2DEG, Z_S_SO2, ZDMS2SO2
REAL(KIND=JPRB)    :: ZDMSMIN
REAL(KIND=JPRB)    :: ZDMSO(KLON), ZGEMU(KLON), ZLATK(KLON)

!-------------------------  Dust constants--------------
INTEGER(KIND=JPIM), PARAMETER :: ISIZE = 11
INTEGER(KIND=JPIM), PARAMETER :: INATS = 12
!INTEGER(KIND=JPIM) , PARAMETER :: INSOIL = 152
INTEGER(KIND=JPIM) , PARAMETER :: INSOIL = 85
INTEGER(KIND=JPIM) , PARAMETER :: INSOILMODE = 3
INTEGER(KIND=JPIM) , PARAMETER :: INBIN = 3
REAL(KIND=JPRB),PARAMETER :: ZZ0SBIS=3.E-05_JPRB   ! en m
REAL(KIND=JPRB), PARAMETER :: ZCSTD = 3.4_JPRB
REAL(KIND=JPRB), PARAMETER :: ZRVONKAR = 0.4_JPRB
REAL(KIND=JPRB), PARAMETER :: ZSIGMAS = 3.0_JPRB
!     Normalization constant
REAL(KIND=JPRB), PARAMETER :: ZRCV=12.62_JPRB
REAL(KIND=JPRB), PARAMETER :: ZLAMBDA=12.0_JPRB
REAL(KIND=JPRB) , DIMENSION(INBIN,2) :: ZTRSIZE
REAL(KIND=JPRB) , DIMENSION(3,12) :: ZMMD , ZPCENT , ZSIGMA
REAL(KIND=JPRB) , DIMENSION(3,12) :: ZMMD2 , ZPCENT2 , ZSIGMA2
REAL(KIND=JPRB), DIMENSION(2,ISIZE)  ::  ZAEROSIZE

!----------------  Other dust variables
REAL(KIND=JPRB) :: ZSTOTAL , ZXK , ZXL , ZXM , ZXN, ZRWI, ZTOTV
REAL(KIND=JPRB) :: ZCOEF, ZDUSTSRC, ZALPHAPROP
REAL(KIND=JPRB) :: ZRE, ZKK, ZDELDP, ZUTH, ZFDP1, ZFDP2, ZDDP, ZTMP, ZZ0M
REAL(KIND=JPRB), DIMENSION(INSOIL) :: ZZ0S
REAL(KIND=JPRB), DIMENSION(INATS,INSOIL) :: ZSREL
REAL(KIND=JPRB), DIMENSION(INSOIL) :: ZSS, ZDP_ARRAY
REAL(KIND=JPRB), DIMENSION(KLON) :: ZCLAY,ZSAND,ZFW,ZFEFF,ZUSTAR, ZSILT
REAL(KIND=JPRB), DIMENSION(INATS,ISIZE,KLON) :: ZRSFROWSUB
REAL(KIND=JPRB), DIMENSION(INATS,KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2)) :: ZRSFROWT
REAL(KIND=JPRB), DIMENSION(KLON,INSOIL) :: ZUSTARTS, ZUSTART,ZFEFFS
REAL(KIND=JPRB), DIMENSION(ISIZE) :: ZFRAC
REAL(KIND=JPRB), DIMENSION(INATS,KLON) :: ZFSOIL, ZFTEX


INTRINSIC ERF

!! REAL(KIND=JPRB)    :: ZDMSSRC(KLON)

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "updcal.intfb.h"
#include "fcttim.func.h"

#include "aer_ssalt.intfb.h"
#include "aer_ssalt_ms.intfb.h"
#include "aer_ssalt_grythe.intfb.h"
#include "aer_ssalt_albert.intfb.h"
#include "abor1.intfb.h"
#include "aer_dmso.intfb.h"
#include "satur.intfb.h"
#include "aer_volce.intfb.h"
#include "aer_stratcl.intfb.h"

IF (LHOOK) CALL DR_HOOK('AER_SRC',0,ZHOOK_HANDLE)


!-----------------------------------------------------------------------
ASSOCIATE(YDCSGLEG=>YDGEOMETRY%YRCSGLEG, YDEAERMAP=>YDMODEL%YRML_PHY_AER%YREAERMAP, YDEAERVOL=>YDMODEL%YRML_PHY_AER%YREAERVOL, &
& YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
& YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
& YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM)

ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, YAERO=>YGFL%YAERO, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, &
& LAERCLIST=>YDEAERATM%LAERCLIST, LAERELVS=>YDEAERATM%LAERELVS, LAERVOL=>YDEAERATM%LAERVOL, NTYPAER=>YDEAERATM%NTYPAER, &
& LAERNITRATE=>YDCOMPO%LAERNITRATE, LAERSOA=>YDCOMPO%LAERSOA,  LAERSOA_COUPLED=>YDCOMPO%LAERSOA_COUPLED, &
& NDUSRCP=>YDEAERMAP%NDUSRCP, RDDUAER=>YDEAERMAP%RDDUAER, RDUSRCP=>YDEAERMAP%RDUSRCP, &
& LOCNDMS=>YDEAERSRC%LOCNDMS, NAERWND=>YDEAERSRC%NAERWND, NDDUST=>YDEAERSRC%NDDUST, NDMSO=>YDEAERSRC%NDMSO, &
& NSSALT=>YDEAERSRC%NSSALT, RAERDUB=>YDEAERSRC%RAERDUB, RCODECA=>YDEAERSRC%RCODECA, RCOVSRA=>YDEAERSRC%RCOVSRA, &
& RDDUSRC=>YDEAERSRC%RDDUSRC, RFCTSS=>YDEAERSRC%RFCTSS, RSIDECA=>YDEAERSRC%RSIDECA,   RSIVSRA=>YDEAERSRC%RSIVSRA, &
& NAERVOLC=>YDEAERVOL%NAERVOLC, NAERVOLE=>YDEAERVOL%NAERVOLE, NVOLDATS=>YDEAERVOL%NVOLDATS, NVOLERUP=>YDEAERVOL%NVOLERUP, &
& NVOLDATE=>YDEAERVOL%NVOLDATE, RAERVOLC=>YDEAERVOL%RAERVOLC, RAERVOLE=>YDEAERVOL%RAERVOLE, &
& LVDFTRAC=>YDEPHY%LVDFTRAC, YSURF=>YDEPHY%YSURF, &
& NLOENG=>YDGEOMETRY%YRGEM%NLOENG, NGLOBALAT=>YDGEOMETRY%YRMP%NGLOBALAT, NSTASS=>YDRIP%NSTASS, RSTATI=>YDRIP%RSTATI, &
& NDRYDEPVEL_DYN=> YDEAERSNK%NDRYDEPVEL_DYN, &
& LAERCHEM=>YGFL%LAERCHEM, &
& RA=>YRCST%RA, RPI=>YRCST%RPI, RDAY=>YRCST%RDAY, RG=>YRCST%RG)


! N.B.: In ECMWF model conventions, flux going upward from the surface 
! are negative
! All surface fluxes PCFLUX in kg m-2 s-1

!-----------------------------------------------------------------------

! N.B.: Security parameters
ZEPSISS=1.E-09_JPRB
ZEPSIDD=1.E-12_JPRB
ZEPSIRA=1.E-06_JPRB
ZEPSISS=0.E+00_JPRB
ZEPSIDD=0.E+00_JPRB
ZEPSIRA=0.E+00_JPRB
ZEPSSNO=1.E-03_JPRB
ZEPSARE=1.E-03_JPRB

Z_S_SO2=0.5_JPRB
ZSO2MSS=64.058E-03_JPRB
ZDMS2SO2=1.0_JPRB
ZDMSMIN = 5.E-11_JPRB
ZRAD2DEG= 180._JPRB/RPI

!-----------------------------------------------------------------------

!*       0.1   TIME AND DATE OF THE MODEL
!              --------------------------
 
IY0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
ID0=NDD(NINDAT)
INC=(NSSSSS + NINT(RSTATI))/NINT(RDAY)
CALL UPDCAL (ID0, IM0, IY0, INC,  IDD, IMM, IYY, IMON, -1)
IMDATE=IYY*10000+IMM*100+IDD

!*       0.2   A LENGTH OF DAY INDEX
!              ---------------------

DO JL=KIDIA,KFDIA
  IGLGLO=NGLOBALAT(KSTGLO+JL-1)
  ZGEMU(JL)=YDCSGLEG%RMU(IGLGLO)                      ! sine of latitude
  ZLAT=ASIN(YDCSGLEG%RMU(IGLGLO))*ZRAD2DEG
  ZLATK(JL)=ZLAT
  ZCOS0=1._JPRB
  ZSIN0=0._JPRB
  PLDAY(JL)=MAX( RSIDECA*ZGEMU(JL)&
   & -RCODECA*RCOVSRA*SQRT(1.0_JPRB-ZGEMU(JL)**2)* ZCOS0&
   & +RCODECA*RSIVSRA*SQRT(1.0_JPRB-ZGEMU(JL)**2)* ZSIN0&
   & ,0.0_JPRB)
  PDMSO(JL)=0._JPRB
  PLISS(JL)=0._JPRB
  PSO2(JL) =0._JPRB
  PTDMS(JL)=0._JPRB         
  PODMS(JL)=0._JPRB
ENDDO

!-----------------------------------------------------------------------

!*       0.5   CONSTANTS AND ACCESS TO DATA
!              ----------------------------

IHTST=20


ZETAH(KIDIA:KFDIA,0)=0._JPRB
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZETA(JL,JK) =PAP(JL,JK) /PAPH(JL,KLEV)
    ZETAH(JL,JK)=PAPH(JL,JK)/PAPH(JL,KLEV)
  ENDDO
ENDDO

ITEST=0
ZDEGRAD= 180._JPRB/RPI
ZDLAT  = 180._JPRB / NDGLG      ! distance in degrees between latitude lines
ZGRDLAT= RPI / NDGLG                          ! distance in radians between latitude lines
ZGRDLAT2=ZGRDLAT*0.55_JPRB

DO JL=KIDIA,KFDIA
  IGLGLO=NGLOBALAT(KSTGLO+JL-1)
  INDLAT(JL)=IGLGLO
  Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
  ZDLON=Z1GP*2.0_JPRB*RPI      ! distance in radians between longitude points on a given latitude line
  ZGRDLON(JL)=ZDLON
  ZLAT=ASIN(YDCSGLEG%RMU(IGLGLO))             ! latitude in radians
  ZA=COS(ZLAT)*SIN(ZDLON/2.0_JPRB)
  ZC=2.0_JPRB * ASIN( MIN(1.0_JPRB,ZA) )
  ZDIST(JL)=RA * ZC /1000.0_JPRB      ! distance in km between longitude points on a given latitude line
 
  ZGLON(JL)=PGELAM(JL)*ZDEGRAD
  ZGLAT(JL)=ZLAT*ZDEGRAD
  ZGDLAT(JL)=ZDLAT
  ZGDLON(JL)=360._JPRB*Z1GP
 ENDDO


ZFAERO (KIDIA:KFDIA,         1:NACTAERO) = 0.0_JPRB
ZAEROK (KIDIA:KFDIA, 1:KLEV, 1:NACTAERO) = PCEN (KIDIA:KFDIA, 1:KLEV, KAERO(1):KAERO(NACTAERO)) 
ZTAEROK(KIDIA:KFDIA, 1:KLEV, 1:NACTAERO) = PTENC(KIDIA:KFDIA, 1:KLEV, KAERO(1):KAERO(NACTAERO))

INBAER=0

!-----------------------------------------------------------------------

!*       0.3   SURFACE WIND VARIABLE RELEVANT FOR SS AND DU EMISSIONS
!              ------------------------------------------------------

IF (NAERWND == 0) THEN
!-- no gust accounted for
  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
ELSEIF (NAERWND == 1) THEN
!-- gust only for SS, 10-m wind for DU
  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
ELSEIF (NAERWND == 2) THEN
!-- gust only for DU, 10-m wind for SS
  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
ELSEIF (NAERWND == 3) THEN
!-- gust for both SS and DU
  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
ENDIF

! correction to account for the decrease of mean wind and gusts with decreasing
! time step
IF (PTSPHY < 1000) THEN
  ! correction not needed for Grythe et al scheme
  IF (NSSALT < 3) THEN
    ZWNDSS(KIDIA:KFDIA)=1.08_JPRB*ZWNDSS(KIDIA:KFDIA)
  ENDIF
ENDIF

!-----------------------------------------------------------------------

!*       1.0   SEA SALT
!              --------

!- Simplistic lifting from surface based on 10-m wind and land-sea mask
ZHSS=MAX(1.0_JPRB,8434._JPRB/1000._JPRB)

IF (NTYPAER(1) /= 0) THEN

!- ECMWF sea salt fluxes come in either 3- or 10-size bins
! 0.03 - 0.50 - 5.0 - 20.
! 0.03 - 0.06 - 0.12 - 0.24 - 0.48 - 0.96 - 1.92 - 3.84 - 7.68 - 15.36 - 30.72

! ZFLX_SSALT is positive in kg m-2 s-1 out of AER_SSALT, 
! but ECMWF conventions have PCFLX as a negative upward flux

  ZFLX_SSALT(KIDIA:KFDIA,:)=0.0_JPRB
  SELECT CASE (NSSALT)

    CASE (1)
      CALL AER_SSALT (YDEAERATM, YDEAERSRC, KIDIA, KFDIA, KLON, NTYPAER(1),&
        & PCI, PLSM, PCLK, ZWNDSS, ZFLX_SSALT)

    CASE (2)
      CALL AER_SSALT_MS (KIDIA, KFDIA, KLON,&
        & PCI, PLSM, PCLK, ZWNDSS, ZFLX_SSALT)

    CASE (3)
      CALL AER_SSALT_GRYTHE (YDEAERATM, YDEAERSNK, KIDIA, KFDIA, KLON,&
        & PCI, PLSM, PCLK, ZWNDSS, PSST, ZRHCL, ZFLX_SSALT)

    CASE (4)
      CALL AER_SSALT_ALBERT (YDEAERATM, YDEAERSRC, KIDIA, KFDIA, KLON,&
        & PCI, PLSM, PCLK, ZWNDSS, PSST, ZRHCL, ZFLX_SSALT)

    CASE DEFAULT
      CALL ABOR1(" SEASALT AEROSOL SCHEME NSSALT WRONG VALUE OR NOT DEFINED")

  END SELECT




  DO JAER=1,NTYPAER(1)
    INBAER=INBAER+1
    DO JL=KIDIA,KFDIA
      ZFLX_SSALT(JL,JAER)=RFCTSS * MAX(ZFLX_SSALT(JL,JAER), ZEPSISS)
    ENDDO
    DO JL=KIDIA,KFDIA
      PCFLX(JL,KAERO(INBAER))=-ZFLX_SSALT(JL,JAER)*1.E+00_JPRB
    ENDDO

!-- if no vertical diffusion, distribute the flux in layers with scale height
!-- between half-levels IHTST-1 and KLEV
    IF (.NOT.LVDFTRAC) THEN
      DO JK=IHTST,KLEV
        DO JL=KIDIA,KFDIA
          ZDETAH(JL,JK) = ZETAH(JL,JK)**ZHSS - ZETAH(JL,JK-1)**ZHSS
          PTENC(JL,JK,KAERO(INBAER)) = PTENC(JL,JK,KAERO(INBAER))+ZFLX_SSALT(JL,JAER)*ZDETAH(JL,JK)
        ENDDO
      ENDDO
    ENDIF
  ENDDO


ENDIF

INBAER=INBAER+NTYPAER(2)

!-----------------------------------------------------------------------

!*       3.0   ORGANIC MATTER (ORGANIC CARBON, POM, SECONDARY ORGANIC MATTER)
!              --------------------------------------------------------------

IF (NTYPAER(3) /= 0) THEN
 INBAER=INBAER+2
ENDIF

!-----------------------------------------------------------------------

!*       4.0   BLACK CARBON
!              ------------

IF (NTYPAER(4) /= 0) THEN
 INBAER=INBAER+2
ENDIF

!-----------------------------------------------------------------------

!*       5.0   SURFACE SULFATE (SO2 --> SO4)
!              -----------------------------

IF (NTYPAER(5) /= 0) THEN
  IF (LAERCHEM) THEN
    INBAER=INBAER+1
  ELSE

    !-- originally, quick fix to produce a flux of oceanic DMS, following Liss & Merlivat, 
    !   1986, in "The Role of Air-Sea Exchange in Geochemical Cycling", 
    !   ed. Buat-Menard, 113-128.

    ! Only add online DMS - prescribed emissions done elsewhere
    IF (NDMSO /= 0 .AND. LOCNDMS) THEN

      CALL AER_DMSO (YDEAERSRC, KIDIA, KFDIA, KLON,&
        & PCI, PDMS, PLDAY, PLSM,  PTS, PWIND,&
        & ZDMSO, PLISS, PTDMS,&
        & PODMS)

      IF (NDMSO == 1) THEN
        PDMSO(KIDIA:KFDIA)=ZDMSO(KIDIA:KFDIA)
      ELSEIF (NDMSO == 2) THEN
        PDMSO(KIDIA:KFDIA)=PODMS(KIDIA:KFDIA)
      ENDIF 
      !-- whatever the source representation (from parametrisation =1 or file =2) 
      !   apply a weighting coefficient representing transfer from DMS to SO2 
      !   30% from Kloster et al., 2005
      PDMSO(KIDIA:KFDIA)=ZDMS2SO2 * PDMSO(KIDIA:KFDIA)

      PCFLX(KIDIA:KFDIA,KAERO(INBAER+2)) = PCFLX(KIDIA:KFDIA,KAERO(INBAER+2)) - PDMSO(KIDIA:KFDIA)
    ENDIF
    INBAER=INBAER+2
  ENDIF
ENDIF


!-----------------------------------------------------------------------

!*       5.5   Secondary Organics Precursors (standalone)
!              -----------------------------

IF (NTYPAER(7) /= 0) THEN
  INBAER = INBAER+3
ENDIF

IF (NTYPAER(8) /= 0) THEN
  IF (LAERSOA_COUPLED) THEN
    INBAER=INBAER+2
  ELSE
    INBAER=INBAER+4
  ENDIF
ENDIF

!-----------------------------------------------------------------------

!*       6.0   FLY ASH AND SO2 FROM VOLCANIC ERUPTIONS
!              ---------------------------------------

!IF (NTYPAER(6) /= 0) THEN
IF (LAERVOL) THEN
! surface flux is set to zero as tendencies are directly updated over the volcano.
  DO JL=KIDIA,KFDIA
    PCFLX(JL,KAERO(INBAER+1))=0._JPRB
  ENDDO

  ZPI2 = 2 * RPI
  ZPIH = RPI/2._JPRB
  ZDLAT  = 180._JPRB / NDGLG
  ZEQUATOR = ZPI2 * RA
  ZLATSQ = ZEQUATOR / (2*NDGLG)  ! length (m) of the latitude side of the grid
  ZGRDLAT= 180._JPRB / NDGLG
  ZGRDLAT2= ZDLAT*0.55_JPRB

  DO JL=KIDIA,KFDIA
    IGLGLO=NGLOBALAT(KSTGLO+JL-1)
!    IGLGLO=MYLATS(JL)
    Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
    ZDLON=Z1GP*2.0_JPRB*RPI
    ZGDLAT(JL)=ZDLAT
    ZGDLON(JL)=360._JPRB*Z1GP 
  ENDDO

!-- for the NAERVOLC continuous volcanoes, no time/height information is assumed  
  IF (NAERVOLC /= 0) THEN
    DO JVOLC=1,NAERVOLC
      ZGELAV=RAERVOLC(JVOLC,1)
      ZGELOV=RAERVOLC(JVOLC,2)
!-- find closest latitude points on the model grid
      ZCIRCLE=ZEQUATOR*COS(ZGELAV*RPI/180._JPRB)

      DO JL=KIDIA,KFDIA
        ZVOLCASH(JL)=0._JPRB
        ZVOLCSO2(JL)=0._JPRB
        IBASPL=KLEV
        ITOPPL=1
        IGVOLC=0

        IGLGLO=NGLOBALAT(KSTGLO+JL-1)
!        IGLGLO=MYLATS(JL)
        Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
        ZGELAT(JL)=ASIN(YDCSGLEG%RMU(IGLGLO))      ! latitude in radians
        ZDLON=Z1GP*2.0_JPRB*RPI                      ! longitude increment in radians
        ZGRDLON(JL)=ZDLON
        ZGRDLON2=0.55_JPRB*ZGRDLON(JL)

        ZGELAA = MIN( 90._JPRB, RAERVOLC(JVOLC,1)+ZGDLAT(JL)*0.55_JPRB)
        ZGELAB = MAX(-90._JPRB, RAERVOLC(JVOLC,1)-ZGDLAT(JL)*0.55_JPRB)
        ZGELOA = MAX(  0._JPRB, RAERVOLC(JVOLC,2)-ZGDLON(JL)*0.55_JPRB)
        ZGELOB = MIN(360._JPRB, RAERVOLC(JVOLC,2)+ZGDLON(JL)*0.55_JPRB)

        IF (ZGLON(JL) >= ZGELOA .AND. ZGLON(JL) < ZGELOB .AND.&
          & ZGLAT(JL) <= ZGELAA .AND. ZGLAT(JL) > ZGELAB ) THEN
          IGVOLC=JL
          ZLONSQ=ZEQUATOR * ZGRDLON(JL) / (2._JPRB*RPI) ! length (m) of the longitude side of the grid 
          ZAREA=ZLATSQ*ZLONSQ
          DO JK=1,KLEV
            IF (RAERVOLC(JVOLC,5) < PALTH(JL,JK) .AND. RAERVOLC(JVOLC,5) >= PALTH(JL,JK+1)) THEN
              IBASPL=JK
            ENDIF
            IF (RAERVOLC(JVOLC,6) < PALTH(JL,JK) .AND. RAERVOLC(JVOLC,6) >= PALTH(JL,JK+1)) THEN
              ITOPPL=JK
            ENDIF
          ENDDO
          IF (IBASPL == ITOPPL) THEN
            ITOPPL=ITOPPL-1
          ENDIF
          ZDEEPLUME=(PAPHI(JL,ITOPPL)-PAPHI(JL,IBASPL))/RG
          ZVOLUME=ZAREA*ZDEEPLUME
          ZVOLCASH(JL)=RAERVOLC(JVOLC,3)*RAERVOLC(JVOLC,7)/ZVOLUME
          ZVOLCSO2(JL)=RAERVOLC(JVOLC,4)*RAERVOLC(JVOLC,8)/ZVOLUME
          DO JK=ITOPPL,IBASPL
!- updating fly ash tendency (should be in variable 13)
            PTENC(JL,JK,KAERO(INBAER+1))=PTENC(JL,JK,KAERO(INBAER+1))+&
              & ZVOLCASH(JL)/PRHO(JL,JK)
!- updating volcanic SO2 tendency (should be in variable 14)
            PTENC(JL,JK,KAERO(INBAER+2))=PTENC(JL,JK,KAERO(INBAER+2))+&
              & ZVOLCSO2(JL)/PRHO(JL,JK)
          ENDDO
        ENDIF
      ENDDO

    ENDDO
  ENDIF
!-- end of contribution by continous volcanoes

  
!-- Explosive volcanoes but no time/height information is assumed  
  IF (NAERVOLE /= 0 .AND. NVOLERUP == 1) THEN    
    DO JVOLE=1,NAERVOLE
      ZGELAV=RAERVOLE(JVOLE,1)
      ZGELOV=RAERVOLE(JVOLE,2)
!-- find closest latitude points on the model grid
      ZCIRCLE=ZEQUATOR*COS(ZGELAV*RPI/180._JPRB)
      IVDATES=NVOLDATS(JVOLE)
      IVDAYS=IVDATES/100
      IVSECNDS=(IVDATES-IVDAYS*100)*3600

      IVDAYE=NVOLDATE(JVOLE)/100
      IVSECNDE=(NVOLDATE(JVOLE)-IVDAYE*100)*3600

      IF ( (IMDATE > IVDAYS .OR. (IMDATE == IVDAYS .AND. NSTASS >= IVSECNDS)) .AND. &
         & (IVDAYE == 0 .OR. IMDATE < IVDAYE .OR. (IMDATE == IVDAYE .AND. NSTASS <  IVSECNDE )) ) THEN

        DO JL=KIDIA,KFDIA
          ZVOLEASH(JL)=0._JPRB
          ZVOLESO2(JL)=0._JPRB
          IBASPL=KLEV
          ITOPPL=1
          IGVOLE=0
          IGLGLO=NGLOBALAT(KSTGLO+JL-1)
!          IGLGLO=MYLATS(JL)
          Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
          ZGELAT(JL)=ASIN(YDCSGLEG%RMU(IGLGLO))      ! latitude in radians
          ZGRDLON(JL)=ZDLON
          ZGRDLON2=1.00_JPRB*ZGRDLON(JL) 
          ZGDLAT(JL)=ZDLAT
          ZGDLON(JL)=360._JPRB*Z1GP 

          ZGELAA = MIN( 90._JPRB, RAERVOLE(JVOLE,1)+ZGDLAT(JL)*0.55_JPRB)
          ZGELAB = MAX(-90._JPRB, RAERVOLE(JVOLE,1)-ZGDLAT(JL)*0.55_JPRB)
          ZGELOA = MAX(  0._JPRB, RAERVOLE(JVOLE,2)-ZGDLON(JL)*0.55_JPRB)
          ZGELOB = MIN(360._JPRB, RAERVOLE(JVOLE,2)+ZGDLON(JL)*0.55_JPRB)

          IF (ZGLON(JL) >= ZGELOA .AND. ZGLON(JL) < ZGELOB .AND.&
            & ZGLAT(JL) <= ZGELAA .AND. ZGLAT(JL) > ZGELAB ) THEN

            IGVOLE=JL
            ZLONSQ = ZEQUATOR * ZGRDLON(JL) / (2._JPRB*RPI)  ! length (m) of the longitude side of the grid 
            ZAREA = ZLATSQ*ZLONSQ
            DO JK=1,KLEV
              IF (RAERVOLE(JVOLE,5) < PALTH(JL,JK) .AND. RAERVOLE(JVOLE,5) >= PALTH(JL,JK+1)) THEN
                IBASPL=JK
              ENDIF
              IF (RAERVOLE(JVOLE,6) < PALTH(JL,JK) .AND. RAERVOLE(JVOLE,6) >= PALTH(JL,JK+1)) THEN
                ITOPPL=JK
              ENDIF
            ENDDO
            IF (IBASPL == ITOPPL) THEN
              ITOPPL=ITOPPL-1
            ENDIF
            ZDEEPLUME=(PAPHI(JL,ITOPPL)-PAPHI(JL,IBASPL))/RG
            ZVOLUME=ZAREA*ZDEEPLUME
            ZVOLEASH(JL)=RAERVOLE(JVOLE,3)*RAERVOLE(JVOLE,7)/ZVOLUME
            ZVOLESO2(JL)=RAERVOLE(JVOLE,4)*RAERVOLE(JVOLE,8)/ZVOLUME

            DO JK=ITOPPL,IBASPL


                PTENC(JL,JK,KAERO(INBAER+1))=PTENC(JL,JK,KAERO(INBAER+1))+&
               & ZVOLEASH(JL)/PRHO(JL,JK) ! kg kg-1 s-1
!- updating volcanic SO2 tendency (should be in variable 15)
              PTENC(JL,JK,KAERO(INBAER+3))=PTENC(JL,JK,KAERO(INBAER+3))+&
               & ZVOLESO2(JL)/PRHO(JL,JK) ! kg kg-1 s-1
                write (*,*) "AERVOLC",JL,JK,ZVOLESO2(JL)/PRHO(JL,JK),PRHO(JL,JK),ZVOLUME,ZAREA,ZLATSQ,ZLONSQ,ZDEEPLUME

            ENDDO

          ENDIF
        ENDDO
      ENDIF
    ENDDO
!-- end of contribution by explosive volcano without emission details

!-- Explosive volcanoes but some volcanoes without and some with time/height information 
  ELSEIF (NAERVOLE /= 0 .AND. NVOLERUP == 2) THEN


     CALL AER_VOLCE&
     &(YDGEOMETRY, YDEAERATM,YDEAERVOL,YDMODEL%YRML_GCONF, &
     & KIDIA, KFDIA, KLON  , KTDIA, KLEV , KSTART, KSTEP, KSTGLO,&
     &KTRAC, KAERO, INBAER,&
     &PALTH, PAPHI, PCEN  , ZGLAT, ZGLON, PRHO  ,&
     &PCFLX, PTENC&
     &)

  ENDIF

  INBAER=INBAER+1

ENDIF
!-----------------------------------------------------------------------

!*       7.0 / 8.0  PSEUDO-STRATOSPHERIC BACKAGROUND AND VOLCANIC AEROSOLS
!                   ------------------------------------------------------
!                   ******* IN PRINCIPLE SHOULD NOT BE USED ANYMORE ******

!-- preliminary configuration: on the first time-step, climatological 
!   stratospheric background and volcanic aerosol optical thicknesses 
!   at 0.55 um are used as initial conditions; they are translated 
!   into mass mixing ratios then advected by the model

IF (KSTEP > KSTART) LAERCLIST = .FALSE.
!IF (LAERCLIST .AND. KSTEP == KSTART .AND. (NTYPAER(7) /= 0 .OR. NTYPAER(8) /= 0)) THEN
IF (LAERCLIST .AND. KSTEP == KSTART .AND. NTYPAER(8) /= 0) THEN
  IFLAG=2
  CALL SATUR (YRTHF, YRCST, KIDIA, KFDIA, KLON  , KTDIA, KLEV, YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
    & PAP, PT   , ZQSAT, IFLAG )  

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZRHCL(JL,JK)=PQ(JL,JK)/ZQSAT(JL,JK)
    ENDDO
  ENDDO
  DO JK=1,KLEV-1
    DO JL=KIDIA,KFDIA
      ZTH(JL,JK)=(PT(JL,JK)*PAP(JL,JK)&
       & *(PAP(JL,JK+1)-PAPH(JL,JK+1))&
       & +PT(JL,JK+1)*PAP(JL,JK+1)*(PAPH(JL,JK+1)-PAP(JL,JK)))&
       & *(1.0_JPRB/(PAPH(JL,JK+1)*(PAP(JL,JK+1)-PAP(JL,JK))))  
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    ZTH(JL,0)=PT(JL,1)-PAP(JL,1)*(PT(JL,1)-ZTH(JL,1))&
     & /(PAP(JL,1)-PAPH(JL,1)) 
    ZTH(JL,KLEV)=PTS(JL)  
  ENDDO

!-- location of reference wavelength 0.55 um
  IWAVL= 9

   CALL AER_STRATCL&
    &(YDGEOMETRY%YRVAB,YDGEOMETRY%YRCVER,YDGEOMETRY%YRDIMV, &
    & YDMODEL%YRML_PHY_RAD%YREAERD,YDMODEL%YRML_PHY_RAD%YRERAD,YDEAERSNK,YDMODEL%YRML_GCONF%YRRIP%YREAERC_TEGEN%RTAEVO,&
    & KIDIA, KFDIA, KLON  , KLEV , IWAVL,&
    & PAPH , ZGEMU, ZRHCL, ZTH,&
    &  ZAEROCLIS )

!-- filling up relevant arrays depending on required aerosol

!-- stratospheric background aerosols
  IF (NTYPAER(7) /= 0) THEN
    DO JL=KIDIA,KFDIA
      PCFLX(JL,KAERO(INBAER+1))=0._JPRB
    ENDDO
!- atmospheric concentration 
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PCEN(JL,JK,KAERO(INBAER+1)) = PCEN(JL,JK,KAERO(INBAER+1)) + ZAEROCLIS(JL,JK,1)
      ENDDO
    ENDDO
    INBAER=INBAER+1
  ENDIF

!-- volcanic background aerosols  
  IF (NTYPAER(8) /= 0) THEN
    DO JL=KIDIA,KFDIA
      PCFLX(JL,KAERO(INBAER+1))=0._JPRB
    ENDDO
!- atmospheric concentration 
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PCEN(JL,JK,KAERO(INBAER+1)) = PCEN(JL,JK,KAERO(INBAER+1)) + ZAEROCLIS(JL,JK,2)
      ENDDO
    ENDDO
    INBAER=INBAER+1
  ENDIF
ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SRC',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SRC

