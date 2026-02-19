! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SRC &
 &( YDGEOMETRY, YDMODEL,KIDIA, KFDIA, KLON , KTDIA, KLEV, KTILES, KSTART, KSTEP ,KSTGLO,&
 &  KSW  , KTRAC, KAERO, &
 &  PALB , PALBD ,PALUVD, PAPHI , &
 &  PAERDEP, PAERLTS, PAERSCC, PAERGUST, PALTH ,&
 &  PDSF, PDSZ, PDMS,&
 &  PAPH , PAP  , PCI  , PGELAM, PGELAT, PGEMU , PFRTI , PHSDFOR,&
 &  PLSM , PCLK , PQ   , PRHO , PSNS  , PT    , PTS   , PTSPHY, PZ0M, KCHEM,&
 &  PWIND, PWS1 , PSST, &
 &  PDMSO, PLDAY, PLISS, PSO2  , PTDMS, PAERFLX, PAERLIF, PAERMAP , PCFLX , PCEN  , PTENC,&
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
!-----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMCST       , ONLY : RA, RPI, RDAY, RG
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
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PALB(KLON), PALBD(KLON,KSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV), PALTH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERDEP(KLON), PAERLTS(KLON), PAERSCC(KLON), PALUVD(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PAERGUST(KLON), PHSDFOR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMS(KLON), PDSF(KLON), PDSZ(KLON)
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAERMAP(KLON,5)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAERFLX(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAERLIF(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
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
REAL(KIND=JPRB) :: ZRWPWP , ZRWSAT , ZSO2MSS
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

!-- map 
INTEGER(KIND=JPIM) :: IFF, ITYPDU
REAL(KIND=JPRB)    :: ZBNDA, ZBNDB, ZBNDC, ZBNDD, ZBNDE, ZBNDF, ZBNDG, ZBNDH, ZBNDI
REAL(KIND=JPRB)    :: ZBNDJ, ZBNDK, ZBNDL, ZBNDM
REAL(KIND=JPRB)    :: ZLAT, ZLON, ZLONE, ZLONGB, ZLONW
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

#include "surf_inq.h"

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
!  Constants for PNabat's scheme

         
!
data ZMMD2/ 690.0_JPRB , 210.0_JPRB , 10.0_JPRB , &
         & 690.0_JPRB , 210.0_JPRB ,  0.0_JPRB , &
         & 690.0_JPRB , 210.0_JPRB ,  0.0_JPRB , &
         & 520.0_JPRB , 100.0_JPRB ,  5.0_JPRB , &
         & 520.0_JPRB ,  75.0_JPRB ,  2.5_JPRB , &
         & 520.0_JPRB ,  75.0_JPRB ,  2.5_JPRB , &
         & 210.0_JPRB ,  75.0_JPRB ,  2.5_JPRB , &
         & 210.0_JPRB ,  50.0_JPRB ,  2.5_JPRB , &
         & 125.0_JPRB ,  50.0_JPRB ,  1.0_JPRB , &
         & 100.0_JPRB ,  10.0_JPRB ,  1.0_JPRB , &
         & 100.0_JPRB ,  10.0_JPRB ,  0.5_JPRB , &
         & 100.0_JPRB ,  10.0_JPRB ,  0.5_JPRB/
! 
data ZSIGMA2/1.6_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &  1.6_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &  1.6_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &  1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &  1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &  1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &  1.7_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &  1.7_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &  1.7_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &  1.8_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &  1.8_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &  1.8_JPRB , 1.8_JPRB , 1.8_JPRB/
! 
data ZPCENT2/1.0_JPRB , 0.00_JPRB , 0.00_JPRB , &
         & 0.90_JPRB , 0.10_JPRB , 0.00_JPRB , &
         & 0.80_JPRB , 0.20_JPRB , 0.00_JPRB , &
         & 0.50_JPRB , 0.35_JPRB , 0.15_JPRB , &
         & 0.45_JPRB , 0.40_JPRB , 0.15_JPRB , & 
         & 0.35_JPRB , 0.50_JPRB , 0.15_JPRB , &
         & 0.30_JPRB , 0.50_JPRB , 0.20_JPRB , & 
         & 0.30_JPRB , 0.50_JPRB , 0.20_JPRB , &
         & 0.20_JPRB , 0.50_JPRB , 0.30_JPRB , &
         & 0.65_JPRB , 0.00_JPRB , 0.35_JPRB , & 
         & 0.60_JPRB , 0.00_JPRB , 0.40_JPRB , &
         & 0.50_JPRB , 0.00_JPRB , 0.50_JPRB/  

!
data ZMMD/1000.0_JPRB , 100.0_JPRB , 10.0_JPRB , &
         &  690.0_JPRB , 100.0_JPRB , 10.0_JPRB , &
         &  520.0_JPRB , 100.0_JPRB ,  5.0_JPRB , &
         &  520.0_JPRB , 100.0_JPRB ,  5.0_JPRB , &
         &  520.0_JPRB ,  75.0_JPRB ,  2.5_JPRB , &
         &  520.0_JPRB ,  75.0_JPRB ,  2.5_JPRB , &
         &  210.0_JPRB ,  75.0_JPRB ,  2.5_JPRB , &
         &  210.0_JPRB ,  50.0_JPRB ,  2.5_JPRB , &
         &  125.0_JPRB ,  50.0_JPRB ,  1.0_JPRB , &
         &  100.0_JPRB ,  10.0_JPRB ,  1.0_JPRB , &
         &  100.0_JPRB ,  10.0_JPRB ,  0.5_JPRB , &
         &  100.0_JPRB ,  10.0_JPRB ,  0.5_JPRB/
! 
data ZSIGMA/1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.6_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.7_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.7_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.7_JPRB , 1.7_JPRB , 1.8_JPRB , &
         &   1.8_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &   1.8_JPRB , 1.8_JPRB , 1.8_JPRB , &
         &   1.8_JPRB , 1.8_JPRB , 1.8_JPRB/
!
data ZPCENT/0.90_JPRB , 0.10_JPRB , 0.00_JPRB , &
         &   0.60_JPRB , 0.30_JPRB , 0.10_JPRB , &
         &   0.60_JPRB , 0.30_JPRB , 0.10_JPRB , &
         &   0.50_JPRB , 0.35_JPRB , 0.15_JPRB , &
         &   0.45_JPRB , 0.40_JPRB , 0.15_JPRB , & 
         &   0.35_JPRB , 0.50_JPRB , 0.15_JPRB , &
         &   0.30_JPRB , 0.50_JPRB , 0.20_JPRB , & 
         &   0.30_JPRB , 0.50_JPRB , 0.20_JPRB , &
         &   0.20_JPRB , 0.50_JPRB , 0.30_JPRB , &
         &   0.65_JPRB , 0.00_JPRB , 0.35_JPRB , & 
         &   0.60_JPRB , 0.00_JPRB , 0.40_JPRB , &
         &   0.50_JPRB , 0.00_JPRB , 0.50_JPRB/ 

!-----------------------------------------------------------------------
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDEAERMAP=>YDMODEL%YRML_PHY_AER%YREAERMAP,YDEAERVOL=>YDMODEL%YRML_PHY_AER%YREAERVOL, &
 & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM)

ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, YAERO=>YGFL%YAERO, &
 & NDGLG=>YDGEOMETRY%YRDIM%NDGLG, &
 & LAERCLIST=>YDEAERATM%LAERCLIST, LAERELVS=>YDEAERATM%LAERELVS, &
 & LAERVOL=>YDEAERATM%LAERVOL, &
 & NDUSRCP=>YDEAERMAP%NDUSRCP, RDDUAER=>YDEAERMAP%RDDUAER, &
 & RDUSRCP=>YDEAERMAP%RDUSRCP, &
 & RRHO_DD=>YDEAERSNK%RRHO_DD, &
 & LOCNDMS=>YDEAERSRC%LOCNDMS, NAERWND=>YDEAERSRC%NAERWND, &
 & LAERNITRATE=>YDCOMPO%LAERNITRATE, &
 & LAERSOA=>YDCOMPO%LAERSOA, &
 & LAERSOA_COUPLED=>YDCOMPO%LAERSOA_COUPLED, &
 & NDDUST=>YDEAERSRC%NDDUST, NDMSO=>YDEAERSRC%NDMSO, NSSALT=>YDEAERSRC%NSSALT, &
 & NTYPAER=>YDEAERATM%NTYPAER, RAERDUB=>YDEAERSRC%RAERDUB, &
 & RCODECA=>YDEAERSRC%RCODECA, RCOVSRA=>YDEAERSRC%RCOVSRA, &
 & RDDUSRC=>YDEAERSRC%RDDUSRC, RDMSMIN=>YDEAERSRC%RDMSMIN, &
 & RFCTSS=>YDEAERSRC%RFCTSS, RSIDECA=>YDEAERSRC%RSIDECA, &
 & RSIVSRA=>YDEAERSRC%RSIVSRA, &
 & NAERVOLC=>YDEAERVOL%NAERVOLC, NAERVOLE=>YDEAERVOL%NAERVOLE, &
 & NVOLDATS=>YDEAERVOL%NVOLDATS, NVOLERUP=>YDEAERVOL%NVOLERUP, &
 & NVOLDATE=>YDEAERVOL%NVOLDATE, &
 & RAERVOLC=>YDEAERVOL%RAERVOLC, RAERVOLE=>YDEAERVOL%RAERVOLE, &
 & NALBEDOSCHEME=>YDEPHY%NALBEDOSCHEME, LVDFTRAC=>YDEPHY%LVDFTRAC, YSURF=>YDEPHY%YSURF, &
 & NLOENG=>YDGEOMETRY%YRGEM%NLOENG, &
 & NGLOBALAT=>YDGEOMETRY%YRMP%NGLOBALAT, &
 & NSTASS=>YDRIP%NSTASS, RHGMT=>YDRIP%RHGMT, RSTATI=>YDRIP%RSTATI,&
 & NDRYDEPVEL_DYN=> YDEAERSNK%NDRYDEPVEL_DYN, &
 & NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, &
 & LAERDUST_NEWBIN=>YDEAERATM%LAERDUST_NEWBIN, &
 & LAERDUSTSIZEVAR=>YDEAERATM%LAERDUSTSIZEVAR, &
 & LAERDUSTSOURCE=>YDEAERATM%LAERDUSTSOURCE, &
 & LAERCHEM=>YGFL%LAERCHEM)




IF (LAERDUST_NEWBIN) THEN
  ZTRSIZE(:,1)=(/0.01_JPRB , 1.0_JPRB , 2.5_JPRB/)
  ZTRSIZE(:,2)=(/1.0_JPRB , 2.5_JPRB ,20.0_JPRB /)
  ZAEROSIZE(1,:)=(/1.0E-08_JPRB , 2.0E-08_JPRB , 4.0E-08_JPRB ,   &
         & 8.0E-08_JPRB , 1.6E-07_JPRB ,3.2E-07_JPRB ,  &
         & 6.4E-07_JPRB , 1.28E-06_JPRB ,         &
         & 2.56E-06_JPRB , 5.12E-06_JPRB , 10.24E-06_JPRB /)
  ZAEROSIZE(2,:)=(/2.0E-08_JPRB , 4.0E-08_JPRB ,  &
         & 8.0E-08_JPRB , 1.6E-07_JPRB ,3.2E-07_JPRB ,  &
         & 6.4E-07_JPRB , 1.28E-06_JPRB ,       &
         & 2.56E-06_JPRB , 5.12E-06_JPRB ,10.24E-06_JPRB ,      &
         & 20.48E-06_JPRB/)
ELSE
  ZTRSIZE(:,1)=(/0.03_JPRB , 0.55_JPRB , 0.9_JPRB/)
  ZTRSIZE(:,2)=(/0.55_JPRB , 0.9_JPRB,20.0_JPRB /)
  ZAEROSIZE(1,:)=(/2.0E-08_JPRB , 4.0E-08_JPRB ,   &
         & 8.0E-08_JPRB , 1.6E-07_JPRB,3.2E-07_JPRB ,  &
         & 5.0E-07_JPRB ,7.0E-07_JPRB,1.0E-6_JPRB, &
         & 2.56E-06_JPRB , 5.12E-06_JPRB , 10.24E-06_JPRB/)
  ZAEROSIZE(2,:)=(/4.0E-08_JPRB ,   &
         & 8.0E-08_JPRB , 1.6E-07_JPRB,3.2E-07_JPRB ,  &
         & 5.0E-07_JPRB ,7.0E-07_JPRB,1.0E-6_JPRB, &
         & 2.56E-06_JPRB , 5.12E-06_JPRB , 10.24E-06_JPRB,20.48E-06_JPRB/)
ENDIF


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

! Initialise variables
ZUSTAR(:)=0.0_JPRB
ZUSTART(:,:)=0.0_JPRB

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

CALL SURF_INQ(YSURF,PRWPWP=ZRWPWP)
CALL SURF_INQ(YSURF,PRWSAT=ZRWSAT)

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


!*       0.6   EMPIRICAL EFFICIENCY FACTORS FOR SOURCES
!              ----------------------------------------

PAERMAP(KIDIA:KFDIA,:) = 0._JPRB
IF (NTYPAER(2) > 0) THEN
 DO JL=KIDIA,KFDIA
  ZLAT=ZGLAT(JL)
  ZLON=ZGLON(JL) 
  ZBNDA= 30._JPRB+(36._JPRB -ZLAT)*14._JPRB/24._JPRB
  ZBNDB= 30._JPRB+(36._JPRB -ZLAT)*40._JPRB/16._JPRB
  ZBNDC= 38._JPRB+(ZLON-124._JPRB)*12._JPRB/29._JPRB
  ZBNDD= 32._JPRB-(ZLON-243._JPRB)* 6._JPRB/21._JPRB

!-- Eastern border Canada/USA
  ZBNDE= 49._JPRB
  IF (ZLON > 268._JPRB .AND. ZLON < 277._JPRB) THEN
    ZBNDE= 49._JPRB-(ZLON-268._JPRB)*7._JPRB/9._JPRB
  ELSEIF (ZLON >= 277._JPRB .AND. ZLON < 285._JPRB) THEN
    ZBNDE= 42._JPRB+(ZLON-277._JPRB)*2._JPRB/8._JPRB
  ELSEIF (ZLON >= 285._JPRB .AND. ZLON < 310._JPRB) THEN
    ZBNDE= 44._JPRB+(ZLON-285._JPRB)*3._JPRB/25._JPRB
  ENDIF

!-- limits Britain
  ZLONGB=-9999._JPRB
  IF (ZLON > 354._JPRB .AND. ZLON < 360._JPRB) THEN
    ZLONGB=ZLON
  ELSEIF (ZLON >= 0._JPRB .AND. ZLON < 3._JPRB) THEN
    ZLONGB=ZLON+360._JPRB
  ENDIF
  ZBNDF= 47._JPRB+(ZLONGB-349._JPRB)*4.5_JPRB/14._JPRB

!-- limits Ireland
  ZBNDG= 61._JPRB-(ZLON-349._JPRB)*7._JPRB/6._JPRB
  ZBNDH= 45._JPRB+(ZLON-349._JPRB)*9._JPRB/6._JPRB

!-- Western border Brazil
  IF (ZLAT <= 4._JPRB .AND. ZLAT > 2._JPRB) THEN
    ZBNDI= 296._JPRB
  ELSEIF (ZLAT <= 2._JPRB .AND. ZLAT > -4._JPRB) THEN
    ZBNDI= 290._JPRB
  ELSEIF (ZLAT <= -4._JPRB .AND. ZLAT > -7._JPRB) THEN
    ZBNDI= 290._JPRB-(-4._JPRB-ZLAT)*4._JPRB/3._JPRB
  ELSEIF (ZLAT <= -7._JPRB .AND. ZLAT > -11._JPRB) THEN
    ZBNDI= 286._JPRB+(-7._JPRB-ZLAT)*4._JPRB/4._JPRB
  ENDIF

  IF (ZLAT <= -11._JPRB .AND. ZLAT > -18._JPRB) THEN
    ZBNDJ= 294._JPRB+(-11._JPRB-ZLAT)*8._JPRB/7._JPRB
  ELSEIF (ZLAT <= -18._JPRB .AND. ZLAT > -27._JPRB) THEN
    ZBNDJ= 302._JPRB+(-18._JPRB-ZLAT)*4._JPRB/9._JPRB
  ELSEIF (ZLAT <= -27._JPRB .AND. ZLAT > -30._JPRB) THEN
    ZBNDJ= 306._JPRB-(-27._JPRB-ZLAT)*3._JPRB/3._JPRB
  ELSEIF (ZLAT <= -30._JPRB .AND. ZLAT >= -34._JPRB) THEN
    ZBNDJ= 303._JPRB+(-30._JPRB-ZLAT)*4._JPRB/4._JPRB
  ENDIF

!-- Northern border India
  IF (ZLON > 70._JPRB .AND. ZLON <= 90._JPRB) THEN
    ZBNDK= 35._JPRB-(ZLON-70._JPRB)*0.5_JPRB
  ENDIF

!-- South border of Asian deserts
  IF (ZLON > 90._JPRB .AND. ZLON <= 135._JPRB) THEN
    ZBNDL= 25._JPRB+(ZLON-90._JPRB)*15._JPRB/45._JPRB
  ENDIF

!-- North limit of the Argentinian pampas
  IF (ZLON > 285._JPRB .AND. ZLON <= 297._JPRB) THEN
    ZBNDM= -42._JPRB+(ZLON-285._JPRB)*6._JPRB/12._JPRB 
  ENDIF

  IFF=0
  ITYPDU=0
 
!-- North America
!  ITYPDU=1
!----- Canada
  IF ( ZLAT >= ZBNDE .AND.&
    &      (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=1
!----- USA
  ELSEIF ( (ZLAT >= ZBNDD .AND. ZLAT < ZBNDE )&
    & .AND. (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=3
  ENDIF
!-- Alaska
  IF ( (ZLAT < 72._JPRB .AND. ZLAT > 52._JPRB)&
    & .AND. (ZLON > 190._JPRB .AND. ZLON <= 219._JPRB) ) THEN
    IFF=2
  ENDIF

!-- Central America
  IF (ZLAT < ZBNDD .AND.&
    &     (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=4
  ENDIF

!-- South America
  IF ( ZLAT < 12._JPRB .AND.&
    &      (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=5
  ENDIF
!-- Brazil
  IF ( (ZLAT <= 4._JPRB .AND. ZLAT > 2._JPRB)&
    &      .AND. (ZLON >= 296._JPRB .AND. ZLON <= 300._JPRB) ) THEN
    IFF=6 
  ENDIF
  IF (ZLAT <= 2._JPRB .AND. ZLAT > -11._JPRB) THEN
      IF (ZLON >= ZBNDI .AND. ZLON < 330._JPRB) THEN
        IFF=6
      ENDIF
  ENDIF
  IF (ZLAT <= -11._JPRB .AND. ZLAT >= -34._JPRB) THEN
      IF (ZLON >= ZBNDJ .AND. ZLON < 330._JPRB) THEN
        IFF=6
      ENDIF
  ENDIF

!-- Western Europe
  IF ( ZLAT > 36._JPRB .AND. ( ZLON >= 330._JPRB .OR. ZLON <= 30._JPRB) ) THEN
    IFF=10
  ENDIF

!----- Iceland
  IF ( (ZLAT < 67._JPRB .AND. ZLAT > 63._JPRB)&
    &      .AND. ( ZLON > 335._JPRB .AND. ZLON < 353._JPRB) ) THEN
    IFF=7
  ENDIF
!----- Britain  
  IF ( (ZLAT < 63._JPRB .AND. ZLAT > ZBNDF)&
    &      .AND. ( ZLON > 354._JPRB .OR. ZLON < 3._JPRB) ) THEN
    IFF=9
  ENDIF
!----- Ireland
  IF ( (ZLAT < ZBNDG .AND. ZLAT > ZBNDH)&
    &      .AND. ( ZLON > 349._JPRB .AND. ZLON < 355._JPRB) ) THEN
    IFF=8
  ENDIF
  ITYPDU=1
  IF ( (IFF >= 1 .AND. IFF <= 10) .OR. IFF == 16) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF        

! ITYPDU=2 
!-- Russia to Urals
  IF ( ZLON > 30._JPRB .AND. ZLON <= 70._JPRB ) THEN
    IF ( ZLAT > 51._JPRB ) THEN
      IFF=11
    ELSEIF ( ZLAT > 36._JPRB ) THEN
      IFF=12
    ENDIF
  ENDIF
  ITYPDU=2
  IF ( IFF >= 11 .AND. IFF <= 12 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF        


!-- Northern Sahara
!  if ( ( zlat <= 36._JPRB .and. zlat >= 21._JPRB) &
!    & .and. ( zlon >= 330._JPRB .or. zlon <= zbnda) ) then
!    iff=13
!-- Northern Sahara (West)
  IF ( ( ZLAT <= 36._JPRB .AND. ZLAT >= 21._JPRB) &
    & .AND. ( ZLON >= 330._JPRB .OR. ZLON < 2._JPRB) ) THEN
    IFF=36
!-- Northern Sahara (East)
  ELSEIF ( ( ZLAT <= 36._JPRB .AND. ZLAT >= 21._JPRB) &
    & .AND. ( ZLON >= 2._JPRB .AND. ZLON <= ZBNDA) ) THEN
    IFF=37
!-- Southern Sahara (West)
  ELSEIF ( ( ZLAT < 21._JPRB .AND. ZLAT >= 13._JPRB) &
    & .AND. ( ZLON >= 350._JPRB .OR. ZLON < 8._JPRB) ) THEN
    IFF=34
!-- Southern Sahara (West2)
  ELSEIF ( ( ZLAT < 21._JPRB .AND. ZLAT >= 13._JPRB) &
    & .AND. ( ZLON >= 330._JPRB .AND. ZLON < 350._JPRB) ) THEN
    IFF=39
!-- Southern Sahara (East)
  ELSEIF ( ( ZLAT < 21._JPRB .AND. ZLAT >= 13._JPRB) &
    & .AND. ( ZLON >= 8._JPRB .AND. ZLON <= ZBNDA) ) THEN
    IFF=35
!-- Central Africa West
  ELSEIF ( (ZLAT < 13._JPRB .AND. ZLAT >= -12._JPRB) &
    & .AND. ( ZLON >= 330._JPRB .OR. ZLON < 20._JPRB) ) THEN
    IFF=38
!-- Central Africa East
  ELSEIF ( (ZLAT < 13._JPRB .AND. ZLAT >= -12._JPRB) &
    & .AND. ( ZLON >= 20._JPRB .AND. ZLON <= 60._JPRB) ) THEN
    IFF=14
!-- Southern Africa
  ELSEIF ( ZLAT < -12._JPRB .AND. ZLAT >= -60._JPRB &
    & .AND. ( ZLON >= 330._JPRB .OR. ZLON <= 60._JPRB) ) THEN
    IFF=15
  ENDIF

  ITYPDU=3
  IF ( (IFF >= 13 .AND. IFF <= 15) .OR. (IFF >= 34 .AND. IFF <= 39) ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF        

!  ITYPDU=4
!-- Australasia
  IF (ZLON > 70._JPRB .AND. ZLON <= 190._JPRB) THEN
    IFF=26

!-- Siberia
    IF (ZLAT <= 90._JPRB .AND. ZLAT > 51._JPRB) THEN
      IFF=16

!-- South Australasia
!---- Tropical Pacific Islands
    ELSEIF ( ZLAT > -10.5_JPRB) THEN
      IFF=27
!---- Australia
    ELSEIF ( ZLAT <= -10.5_JPRB .AND. ZLAT >= -60._JPRB) THEN
      IFF=28
    ENDIF
  ENDIF
  ITYPDU=4
  IF ( IFF >= 26 .AND. IFF <= 28 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF    

! ITYPDU=5
!-- Asian deserts
  IF ( ZLON > 90._JPRB .AND. ZLON <= 135._JPRB ) THEN
    IF ( ZLAT <= 51._JPRB .AND. ZLAT > ZBNDL ) THEN
      IFF=17
    ENDIF
  ENDIF

!-- Saudi Arabia
  IF (ZLAT <= 36._JPRB .AND. ZLAT >= 12._JPRB) THEN
    IF ( ZLON > ZBNDA .AND. ZLON < ZBNDB) THEN
      IFF=18
    ENDIF
  ENDIF
!-- Irak, Iran, Pakistan
  IF (ZLAT <= 36._JPRB   .AND. ZLAT >= 20._JPRB) THEN
    IF ( ZLON > ZBNDB .AND. ZLON < 70._JPRB) THEN
      IFF=19
    ENDIF
  ENDIF

!-- Central Asia and India
  IF ( ZLON > 70._JPRB .AND. ZLON <= 90._JPRB) THEN
!----- Central Asia: Taklamakan
    IF (ZLAT <= 43._JPRB .AND. ZLAT >= ZBNDK) THEN
      IFF=20
!----- India
    ELSEIF (ZLAT <= ZBNDK .AND. ZLAT > 7._JPRB) THEN
      IFF=21
    ENDIF
  ENDIF
!-- other Gobi(s) in South Mongolia and Central China
  IF ( ZLAT <= 49._JPRB .AND. ZLAT > 35._JPRB) THEN
    IF (ZLON > 90._JPRB .AND. ZLON <= 110._JPRB)  THEN
      IFF=22
    ELSEIF (ZLON > 110._JPRB .AND. ZLON <= 125._JPRB) THEN
      IFF=23
    ENDIF
  ENDIF

!-- South China  
  IF ( ZLON > 90._JPRB .AND. ZLON <= 135._JPRB ) THEN
    IF ( ZLAT <= ZBNDL .AND. ZLAT > 7._JPRB ) THEN
      IFF=24
    ENDIF
  ENDIF
  ITYPDU=5
  IF ( IFF >= 17 .AND. IFF <= 24 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF    


! ITYPDU=7
!-- Japan and S.Korea
  IF (ZLON > 124._JPRB .AND. ZLON < 153._JPRB) THEN
    IF ( ZLAT > 24._JPRB .AND. ZLAT < ZBNDC ) THEN
      IFF=25
    ENDIF
  ENDIF

!-- Greenland
  IF (ZLAT > 50._JPRB) THEN
    ZINCLAT=(90._JPRB-ZLAT)/40._JPRB*45._JPRB
    ZLONW=270._JPRB +ZINCLAT
    ZLONE=360._JPRB -ZINCLAT
    IF ( ZLON > ZLONW .AND. ZLON <  ZLONE ) THEN
      IFF=29
    ENDIF
  ENDIF

!-- Antarctica
  IF (ZLAT < -60._JPRB) THEN
    IFF=30
  ENDIF
  ITYPDU=7
  IF ( (IFF >= 29 .AND. IFF <= 30) .OR. IFF == 25 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF

!-- awaiting a proper recoding, new areas are set between iff=31 and 35
! ITYPDU=6 
  IF ( ZLON > 285._JPRB .AND. ZLON < 295._JPRB) THEN
!- Atacama desert and Salar de Uyuni
    IF ( ZLAT < -16._JPRB  .AND. ZLAT > -28._JPRB) THEN
      IFF=31
    ENDIF
!- Salar de Pipanaco and other small ones
    IF ( ZLAT <= -28._JPRB ) THEN 
      IF ( ZLAT > ZBNDM) THEN
        IFF=32
      ENDIF
    ENDIF
  ENDIF
!- Argentinian pampas
  IF ( ZLON > 285._JPRB .AND. ZLON < 297._JPRB ) THEN
    IF ( ZLAT <= ZBNDM ) THEN
      IFF=33
    ENDIF
  ENDIF
  ITYPDU=6
  IF ( IFF >= 31 .AND. IFF <= 33 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF
       
  IF (IFF /= 0) THEN
    PAERMAP(JL,1)=IFF*PLSM(JL)                                   ! area index - diagnostic only
    PAERMAP(JL,3)=RDUSRCP(NDUSRCP(IFF),1)                        ! reference speed - ACTUALLY USED BELOW for Ginoux
    PAERMAP(JL,4)=RDUSRCP(NDUSRCP(IFF),2)                        ! reference particule radius - ACTUALLY USED BELOW for Ginoux
    ZDUEMPOT(JL,:)=RDDUAER(IFF)  ! dust emission potential factor
  ELSE
    WRITE(NULOUT,FMT='(''aer_src: Unassigned grid for Lat,Lon='',2F8.2)') ZLAT,ZLON
    PAERMAP(JL,1)=0._JPRB
    PAERMAP(JL,3)=6._JPRB ! reference speed for unknown point
    PAERMAP(JL,4)=5._JPRB ! reference radius for unknown point
    ZDUEMPOT(JL,:)=1._JPRB
  ENDIF
  IF (NDDUST /= 2 .AND. NDDUST /= 4) THEN
    ! (including land-sea mask for Ginoux, not for Nabat)
    ZDUEMPOT(JL,:) = ZDUEMPOT(JL,:) * RDDUSRC(1:NTYPAER(2))*PLSM(JL)
  ENDIF
  PAERMAP(JL,2)=ZDUEMPOT(JL,1)                                 ! for diagnostics only
 ENDDO
ENDIF


ZFAERO (KIDIA:KFDIA,         1:NACTAERO) = 0.0_JPRB
ZAEROK (KIDIA:KFDIA, 1:KLEV, 1:NACTAERO) = PCEN (KIDIA:KFDIA, 1:KLEV, KAERO(1):KAERO(NACTAERO)) 
ZTAEROK(KIDIA:KFDIA, 1:KLEV, 1:NACTAERO) = PTENC(KIDIA:KFDIA, 1:KLEV, KAERO(1):KAERO(NACTAERO))

INBAER=0

!-----------------------------------------------------------------------

!*       0.3   SURFACE WIND VARIABLE RELEVANT FOR SS AND DU EMISSIONS
!              ------------------------------------------------------

IF (NAERWND == 0) THEN
!-- no gust accounted for
  ZWNDDU(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
ELSEIF (NAERWND == 1) THEN
!-- gust only for SS, 10-m wind for DU
  ZWNDDU(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
ELSEIF (NAERWND == 2) THEN
!-- gust only for DU, 10-m wind for SS
  ZWNDDU(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
ELSEIF (NAERWND == 3) THEN
!-- gust for both SS and DU
  ZWNDDU(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
ENDIF

! correction to account for the decrease of mean wind and gusts with decreasing
! time step
IF (PTSPHY < 1000) THEN
  IF (NDDUST == 3) THEN
    ZWNDDU(KIDIA:KFDIA)=1.06_JPRB*ZWNDDU(KIDIA:KFDIA)
  ENDIF
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

!-----------------------------------------------------------------------

!*       2.0   DESERT DUST
!              -----------

!- Simplistic lifting from surface based on 10-m wind and surface albedo
ZHDD=MAX(1.0_JPRB,8434._JPRB/1000._JPRB)

PAERLIF(KIDIA:KFDIA,:)=0._JPRB
PAERFLX(KIDIA:KFDIA,:)=0._JPRB
ZFLX_SDUST(KIDIA:KFDIA,:)=0._JPRB
ZFSOIL(:,:)=0.0_JPRB


IF (NTYPAER(2) /= 0) THEN
 !-----------------------------------------------

 IF (NDDUST == 2 .OR. NDDUST == 4 ) THEN ! case P. Nabat formulation

 !-----------------------------------------------------------------------------------------
 ! -------------- Definition of soil textures
 !
 !  CLAY : CLAY >= 0.40 SILT < 0.40 SAND < 0.45
 ! SANDY CLAY : CLAY >= 0.36 SAND >= 0.45 
 ! SILTY CLAY : CLAY >= 0.40 SILT >= 0.40 
 ! SILT : SILT >= 0.8 CLAY < 0.12
 ! SAND : SAND >= 0.3*CLAY + 0.87
 ! SANDY CLAY LOAM : CLAY >= 0.28  CLAY < 0.36 SAND >= 0.45 | CLAY >= 0.20 CLAY <
 ! 0.28 SILT < 0.28
 ! SILTY CLAY LOAM : CLAY >= 0.28 CLAY < 0.40 SAND < 0.20
 ! CLAY LOAM : CLAY >= 0.28 CLAY < 0.40 SAND >= 0.20 SAND < 0.45
 ! SILT LOAM : SILT >= 0.8 CLAY >= 0.12 | SILT >= 0.5 SILT < 0.8 CLAY < 0.28
 ! LOAMY SAND : SAND >= CLAY + 0.7 SAND < 0.3*CLAY + 0.87
 ! SANDY LOAM : SAND >= 0.52  CLAY < 0.20 | SAND >= (0.5 - CLAY)  CLAY < 0.07
 ! LOAM : CLAY >= 0.20 CLAY < 0.28 SILT >= 0.28 SILT < 0.5 | SAND >= (0.5 - CLAY)
 ! CLAY < 0.20
 ZFTEX(:,:)=0.0

 DO JL=KIDIA,KFDIA
  !
  ZSAND(JL)=PAERLTS(JL)
  ZCLAY(JL)=PAERSCC(JL)
  ZSILT(JL) = 1.0_JPRB - ZCLAY(JL) - ZSAND(JL)

  IF (ZSILT(JL) <= 0._JPRB) THEN
     ZSILT(JL) = 0.0_JPRB
  ENDIF
  IF (ZSILT(JL) + 1.5_JPRB*ZCLAY(JL)<0.15_JPRB) THEN  ! Sand
       ZFTEX(1,JL)=1.0_JPRB
  ELSEIF (ZSILT(JL) + 2._JPRB*ZCLAY(JL)<0.30_JPRB) THEN  ! Loamy sand
       ZFTEX(2,JL)=1.0_JPRB
  ELSEIF ( ( ZCLAY(JL)>=0.07_JPRB .AND. ZCLAY(JL)<0.20_JPRB .AND. &
         &   ZSAND(JL)>0.52_JPRB .AND. ZSILT(JL)+2._JPRB*ZCLAY(JL)>=0.30_JPRB ) &
         & .OR. &
         & ( ZCLAY(JL)<0.07_JPRB .AND. ZSILT(JL)<0.50_JPRB .AND. &
         &   (ZSILT(JL)+2._JPRB*ZCLAY(JL))>=0.30_JPRB) ) THEN  !_JPRB .ANDy Loam
       ZFTEX(3,JL)=1.0_JPRB
  ELSEIF ((ZSILT(JL)>=0.50_JPRB .AND. ZCLAY(JL)>=0.12_JPRB .AND. ZCLAY(JL)<0.27_JPRB) .OR. &
        &       (ZSILT(JL)>=0.50_JPRB .AND. ZSILT(JL)<0.80_JPRB .AND. ZCLAY(JL)<0.12_JPRB)) THEN
! Silt Loam
       ZFTEX(4,JL)=1.0_JPRB
  ELSEIF (ZSILT(JL)>=0.80_JPRB .AND. ZCLAY(JL)<0.12_JPRB) THEN  ! Silt
       ZFTEX(5,JL)=1.0_JPRB
  ELSEIF (ZCLAY(JL)>=0.07_JPRB .AND. ZCLAY(JL)<0.27_JPRB .AND. &
        & ZSILT(JL)>=0.28_JPRB .AND. ZSILT(JL)<0.50_JPRB .AND. &
        & ZSAND(JL)<=0.52_JPRB) THEN
       ZFTEX(6,JL)=1.0_JPRB                             ! Loam
  ELSEIF (ZCLAY(JL)>=0.20_JPRB .AND. ZCLAY(JL)<0.35_JPRB .AND. &
        & ZSILT(JL)<0.28_JPRB .AND. ZSAND(JL)>0.45_JPRB) THEN !_JPRB .ANDy Clay Loam
       ZFTEX(7,JL)=1.0_JPRB
  ELSEIF (ZCLAY(JL)>=0.27_JPRB .AND. ZCLAY(JL)<0.40_JPRB .AND. ZSAND(JL)<=0.20_JPRB) THEN !Silty Clay Loam 
       ZFTEX(8,JL)=1.0_JPRB
  ELSEIF (ZCLAY(JL)>=0.27_JPRB .AND. ZCLAY(JL)<0.40_JPRB .AND. ZSAND(JL)>0.20_JPRB .AND. ZSAND(JL)<=0.45_JPRB) THEN ! Clay Loam
       ZFTEX(9,JL)=1.0_JPRB
  ELSEIF (ZCLAY(JL)>=0.35_JPRB .AND. ZSAND(JL)>0.45_JPRB) THEN  !Sandy Clay
       ZFTEX(10,JL)=1.0_JPRB
  ELSEIF (ZCLAY(JL)>=0.40_JPRB .AND. ZSILT(JL)>=0.40_JPRB) THEN  ! Silty Clay
       ZFTEX(11,JL)=1.0_JPRB
  ELSEIF (ZCLAY(JL)>=0.40_JPRB .AND. ZSAND(JL)<=0.45_JPRB .AND. ZSILT(JL)<0.40_JPRB) THEN  !Clay
       ZFTEX(12,JL)=1.0_JPRB
  ENDIF
 ENDDO

 ! Calculation of ZDP_ARRAY
 ZDELDP = 0.0460517018598807_JPRB        ! zdp_array(2)-zdp_array(1) zss(nsoil)dln(Dp)
 ZDP_ARRAY(1)=0.0001_JPRB  ! in cm
 DO JS=2,INSOIL
      ZDP_ARRAY(JS) = ZDP_ARRAY(JS-1)*exp(ZDELDP)
      !write (*,*) "ZDELDP",JS,ZDP_ARRAY(JS),((EXP(ZDELDP))**(JS-1))*(EXP(ZDELDP)-1)*ZDP_ARRAY(1)
 ENDDO

 ! Determination of ZSREL=dSrel in MB95
 ! equations 29-30-31 in MB95
 ! ZSS=dS(Dp)
 DO JNT = 1 , INATS
  ZSS(:) =0.0
  ZSTOTAL = 0.0
  DO JS = 1 , INSOIL          !soil size segregatoin no
    ZDDP=((EXP(ZDELDP))**(JS-1))*(EXP(ZDELDP)-1)*ZDP_ARRAY(1)
    DO JNM = 1 , INSOILMODE       !soil mode = 3
      IF ( ZPCENT(JNM,JNT)>0 .AND. ZSIGMA(JNM,JNT)>0 ) THEN
        ZXK = ZPCENT(JNM,JNT)/(sqrt(2*RPI)*LOG(ZSIGMA(JNM,JNT)))
        ZXL =((LOG(ZDP_ARRAY(JS))-LOG(ZMMD(JNM,JNT)*1.E-4))**2)/(2.0*(LOG(ZSIGMA(JNM,JNT)))**2)
        ZXM = ZXK*EXP(-ZXL)
      ELSE
        ZXM = 0.0
      ENDIF
        ZXN = 2.0*RPI*(ZDP_ARRAY(JS)*1.E4)**2
        ZSS(JS) = ZSS(JS) + (ZXM*ZDELDP*ZXN)
    ENDDO
    ZSTOTAL = ZSTOTAL + ZSS(JS)*ZDDP
    !write (*,*) "ZSRELCOMP",JNT,JS,ZSS(JS),ZDDP,ZSS(JS)*ZDDP,ZSTOTAL
  ENDDO
  DO JS = 1 , INSOIL
    ZDDP=((EXP(ZDELDP))**(JS-1))*(EXP(ZDELDP)-1)*ZDP_ARRAY(1)
    IF (ZSTOTAL > 0.0 ) THEN
      ZSREL(JNT,JS) = ZSS(JS)/ZSTOTAL
      !write (*,*) "ZSREL",JS,JNT,ZSREL(JNT,JS),ZSS(JS),ZSTOTAL
    ENDIF
  ENDDO
 ENDDO                    ! soil types   

 ! Emission distribution (Kok, 2011)
 ZTOTV = 0.0
 DO JN=1,ISIZE
   ZRWI = (ZAEROSIZE(1,JN)+ZAEROSIZE(2,JN))/2.0*1.E6
   ZFRAC(JN) = ZRWI/ZRCV*(1+ERF(LOG(ZRWI/ZCSTD)/SQRT(2.0)/LOG(ZSIGMAS)))*EXP(-(ZRWI/ZLAMBDA)**3)
   !write (*,*) "KOK",JN,ZAEROSIZE(1,JN),ZAEROSIZE(2,JN),ZRWI,ZFRAC(JN)
 !see Kok (2011)
   ZTOTV = ZTOTV + ZFRAC(JN)
 ENDDO

 DO JN=1,ISIZE
  ZFRAC(JN) = ZFRAC(JN)/ZTOTV
  !if ( ZFRAC(JN).lt.1.E-9 ) ZFRAC(JN) = 0.0
 ENDDO

   !-- PRELIMINARY FOR TESTING

   !-- surface source of dust is assumed if fraction of bare soil is > 0.1
   !and/or UVis albedo > 0.11
   !                                        10m wind > threshold = f(soil
   !                                        wetness, particle radius)
   !                                        no snow, and below the wilting point
   !                                        (0.171)
 DO JL=KIDIA,KFDIA
  LLDUST(JL)=.TRUE.
  ZSAND(JL)=PAERLTS(JL)*100._JPRB
  ZCLAY(JL)=PAERSCC(JL)*100._JPRB

      !-------------------------------------------------------------------
      ! ---  MB95 - Modifications by P.Nabat (10/2013)

      !      IF (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND.
      !      PFRTI(JL,8) >= 0.90_JPRB &
      !      Sensitivity test adopted
!  LLDDACT(JL)= (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND.PFRTI(JL,8) >= 0.70_JPRB &
!  &   .AND. PWS1(JL) <= 2._JPRB * ZRWPWP  .AND. ZWNDDU(JL).gt.0.)
  !.AND. PHSDFOR(JL) <= 50._JPRB 

  LLDDACT(JL)= (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) < ZEPSSNO .AND. PHSDFOR(JL) <= 50._JPRB .AND. PALB(JL) <0.52_JPRB .AND.&
        & PFRTI(JL,2) < ZEPSARE .AND. PFRTI(JL,3) < ZEPSARE .AND.&! no ice, no wet skin
        & PFRTI(JL,5) < ZEPSARE .AND. PFRTI(JL,6) < ZEPSARE .AND.&! no snow under bare soil-low veg, no dry high veg
        & PFRTI(JL,7) < ZEPSARE .AND.&! no snow under high veg
       & PFRTI(JL,8) > 0.8_JPRB .AND. PFRTI(JL,4) < 0.5_JPRB )
  IF (LAERDUSTSOURCE) THEN
    LLDDACT(JL)= (PLSM(JL) >= 0.19_JPRB .AND. PSNS(JL) < ZEPSSNO .AND. &
    & PFRTI(JL,2) < ZEPSARE .AND. PFRTI(JL,3) < ZEPSARE .AND.&! no ice, no wet skin
    & PFRTI(JL,5) < ZEPSARE .AND. PFRTI(JL,6) < 0.3_JPRB .AND. &! no snow under bare soil-low veg, no dry high veg
    & PFRTI(JL,7) < ZEPSARE )
  ENDIF

  IF (LLDDACT(JL)) THEN
    IF (NDDUST == 4) THEN
      DO JS=1,INSOIL
     ! Calculation of Z0s : roughness length of smooth erodible surfaces
     ! Equ (18) in Marticorena et al. (1997)
     ZZ0S(JS) = (ZDP_ARRAY(JS)*1.E-2)/30.
          ! ZFEFFS = correction factor for the effect of surface roughness
          ! Equ (3) in Zakey et al. (2006)
          ! Equ (20) in Marticorena and Bergametti (1995)
          ! Equ (2) in Laurent et al. (2008)
          ! Equ (4) in Menut et al. (2013) with error in units (Z0s should be in
          ! m)
          !ZFEFFS(JL,JS)=MIN(1.0,MAX(0.001,1-(LOG(PGZ0(JL)/(200.*ZZ0S(JS)))/LOG(0.35*(0.1/ZZ0S(JS))**0.8))))
          ZFEFFS(JL,JS)=MIN(1.0_JPRB,MAX(0.001_JPRB,1._JPRB-(LOG(PZ0M(JL)/(100._JPRB*ZZ0S(JS))) &
           &                                                 /LOG(0.35_JPRB*(0.1_JPRB/ZZ0S(JS))**0.8_JPRB))))
      ENDDO
      ! Why do we overwrite this with the value for one particular JS index?
      PAERMAP(JL,2)=ZFEFFS(JL,40)
    ELSE   ! NDDUST==2
      ZFEFF(JL)=MIN(1.0_JPRB,MAX(0.001_JPRB,1._JPRB-(LOG(0.0002_JPRB/ZZ0SBIS)/LOG(0.35_JPRB*(0.1_JPRB/ZZ0SBIS)**0.8_JPRB))))
      PAERMAP(JL,2)=ZFEFF(JL)
    ENDIF
    ! ZFW = factor for the effect of soil moisture = Wetu*T/DryU*T
    ! ZSWETN=w "total soil moisture" in Fecan et al. (1999)
    ! w="prognostic surface volumetric soil moisture calculated by BATS"
    ! in Zakey et al.(2006)
    ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB,(PWS1(JL)-ZRWPWP)/(ZRWSAT-ZRWPWP) ) )
    ! ZSWETN2=w' in Fecan et al. (1999) : empirical equation (14/15)
    ! w'=residual soil moisture (Laurent et al. 2008)
    ZSWETN2 = 0.0014_JPRB*ZCLAY(JL)**2+0.17_JPRB*ZCLAY(JL)   ! 0.0015 in Zakey et al. (2006) ?
    IF (ZSWETN <= ZSWETN2) THEN
      ZFW(JL)=1._JPRB
    ELSE
      ZFW(JL)=(1._JPRB+1.21_JPRB*(ZSWETN-ZSWETN2)**0.68_JPRB)**0.5_JPRB   ! Fecan et al.(1999) equ (15)
    ENDIF
  ENDIF
 ENDDO

   ! Calculation of ZUSTARTS = ideal minimum threshold friction velocity
   !   and ZUSTART = threshold wind friction velocity
   DO JS = 1,INSOIL
     ! Marticorena and Bergametti (1995) :
     ! ZRE= Reynolds number (equation 5)
     ZRE = 1331.647_JPRB*ZDP_ARRAY(JS)**1.561228_JPRB+0.38194_JPRB    ! With ZDP_ARRAY in cm

     DO JL=KIDIA,KFDIA
       IF (LLDDACT(JL)) THEN
          ! Marticorena and Bergametti (1995), ZKK=K in equation (4)
          ZKK = (RRHO_DD(1)/PRHO(JL,KLEV) * RG * ZDP_ARRAY(JS)*1.E-2_JPRB * & ! RRHO_DD and PRHO in kg.m-3
          &        (1._JPRB+0.006_JPRB/(0.001_JPRB * RRHO_DD(1) * 100._JPRB*RG &
          &                              *(ZDP_ARRAY(JS))**2.5_JPRB)))**0.5_JPRB  !Put RG in cm.s-2 ???
          IF (ZRE > 10._JPRB) THEN
             ! MB95 equation (7)
             ZUSTARTS(JL,JS)=0.129_JPRB*ZKK*(1._JPRB-0.0858_JPRB*exp(-0.0617_JPRB*(ZRE-10._JPRB)))
          ELSE
             ! MB95 equation (6)
             ZUSTARTS(JL,JS)=0.129_JPRB*ZKK/((1.928_JPRB*(ZRE**0.092_JPRB)-1._JPRB)**0.5_JPRB)
          ENDIF
          ! MB95 equation (10/17/21)
          ! Zakey et al. (2006) equ (1)
          ! Laurent et al. (2008) equ (1)
          ! Menut et al. (2013) equ (3)
          IF (NDDUST == 4) THEN
            ZUSTART(JL,JS) = ZUSTARTS(JL,JS)/ZFEFFS(JL,JS)*ZFW(JL)
          ELSE
            ZUSTART(JL,JS) = ZUSTARTS(JL,JS)/ZFEFF(JL)*ZFW(JL)
          ENDIF
       ENDIF
     ENDDO
   ENDDO

 ! Saltation flux (Marticorena and Bergametti, 1995)
   IF (LAERDUSTSOURCE) THEN
     ZCOEF=3.E-2_JPRB
   ELSE
     IF (.NOT.LAERDUST_NEWBIN) THEN
       ZCOEF = 1.3E-2_JPRB
       IF (NDRYDEPVEL_DYN > 0) THEN
         ZCOEF = 1.7E-2_JPRB
       ENDIF
       IF (NDDUST == 4) THEN
         ZCOEF = ZCOEF*3._JPRB
       ENDIF
     ELSE
       ZCOEF = 1.7E-2_JPRB
       IF (NDDUST == 4) THEN
         ZCOEF = 4.8E-2_JPRB
         IF (NDRYDEPVEL_DYN >0) THEN
           ZCOEF = 5.4E-2_JPRB
         ENDIF
       ENDIF
     ENDIF
     ZCOEF=ZCOEF*1.2_JPRB
   ENDIF

   
   DO JL=KIDIA,KFDIA
    IF (LLDDACT(JL)) THEN
     ZZ0M=MIN(PZ0M(JL),0.5_JPRB)
     ZUSTAR(JL) = ZWNDDU(JL)*ZRVONKAR/LOG(10._JPRB/ZZ0M)
    ENDIF
   ENDDO

   DO JS=1,INSOIL
     DO JL=KIDIA,KFDIA
      IF (LLDDACT(JL) .AND. ZUSTAR(JL)>ZUSTART(JL,JS)) THEN
       ZUTH = ZUSTART(JL,JS)/ZUSTAR(JL)
       ZFDP1 = ZUSTAR(JL)**3._JPRB * MAX(0._JPRB,(1.0_JPRB-ZUTH**2._JPRB))
       ZFDP2 = (1._JPRB+ZUTH)*PRHO(JL,KLEV)/RG

       ! Conversion from horizontal to vertical flux (sandblasting)
       ! Gilette (1979) - equ (47) in MB95       
       IF (ZCLAY(JL) <= 17.0_JPRB) THEN  !before correction in Eastern Africa: 20.0
           ZALPHAPROP = ZCOEF*10._JPRB**(0.134_JPRB*ZCLAY(JL)-6.0_JPRB)
           !ZALPHAPROP = 10._JPRB**(0.134_JPRB*ZCLAY(JL)-6.0_JPRB)
       ELSE
          !ZALPHAPROP = ZCOEF*10**(-0.1*ZCLAY(JL)-1.2)
          ZALPHAPROP = ZCOEF*10._JPRB**(-0.09_JPRB*ZCLAY(JL)-2.19_JPRB)
          !ZALPHAPROP = 10._JPRB**(-0.09_JPRB*ZCLAY(JL)-2.19_JPRB)
       ENDIF

       IF (ZSAND(JL) > 0._JPRB) THEN
         DO JNT=1,INATS
           ZTMP=ZFSOIL(JNT,JL)
           !ZFSOIL(JNT,JL) = ZFSOIL(JNT,JL) + ZALPHAPROP * MIN(ZSREL(JNT,JS),1.0_JPRB) * ZFDP1* ZFDP2 * &
           ZFSOIL(JNT,JL) = ZFSOIL(JNT,JL) + ZALPHAPROP * ZSREL(JNT,JS) * ZFDP1* ZFDP2 * &
            &                 ((EXP(ZDELDP))**(JS-1))*(EXP(ZDELDP)-1)*ZDP_ARRAY(1)
         IF (ZFSOIL(JNT,JL) > 1.E-8_JPRB) THEN
         !Write (*,*) "AERSRC DUST, JL",JL,JNT,JS,ZCLAY(JL),ZALPHAPROP,ZFSOIL(JNT,JL),ZTMP,ZSREL(JNT,JS),ZFDP1,ZFDP2,&
         ! & ZUSTAR(JL),ZUSTART(JL,JS),(EXP(ZDELDP))**(JS-1)*(EXP(ZDELDP)-1)*ZDP_ARRAY(1),PAERUST(JL)
         ENDIF
         ENDDO
       ENDIF
      ENDIF
     ENDDO
   ENDDO
   !!!!!!!!!!!!!!
              ! Calculate vertical fluxes (Kok distribution)
   ZRSFROWT(:,:,:) = 0.0_JPRB
   DO JAER=1,NTYPAER(2)
     DO JL=KIDIA,KFDIA
      IF (LLDDACT(JL)) THEN
       DO JN=1,ISIZE
         DO JNT=1,INATS
             ZRSFROWSUB(JNT,JN,JL) = ZFSOIL(JNT,JL)*ZFRAC(JN)
             ZRWI = (ZAEROSIZE(1,JN)+ZAEROSIZE(2,JN))/2.0_JPRB*1.E6_JPRB
             ! write (*,*) "FINAL D",JAER,JN,ZRWI,ZTRSIZE(JAER,1),ZTRSIZE(JAER,2),ZFRAC(JN)
             IF (ZRWI >= ZTRSIZE(JAER,1) .AND. ZRWI < ZTRSIZE(JAER,2)) THEN
                        ZRSFROWT(JNT,JL,JAER) = ZRSFROWT(JNT,JL,JAER) + ZRSFROWSUB(JNT,JN,JL)
             ENDIF
          ENDDO
        ENDDO
       ENDIF
     ENDDO

    DO JL=KIDIA,KFDIA
      ZFLX_SDUST(JL,JAER)=0._JPRB
      IF (LLDDACT(JL)) THEN

        ZDUSTSRC=1._JPRB

        DO JNT=1,INATS
             ZFLX_SDUST(JL,JAER) = ZFLX_SDUST(JL,JAER) + &
              & ZRSFROWT(JNT,JL,JAER)*ZFTEX(JNT,JL)*PFRTI(JL,8)*(1._JPRB-PSNS(JL))*PLSM(JL)*ZDUSTSRC

        ENDDO

       ! Use of Ginoux (MODIS deep blue DOD) dust source function
       IF (LAERDUSTSOURCE) THEN
         ZFLX_SDUST(JL,JAER)=ZFLX_SDUST(JL,JAER)*PDSF(JL)
       ELSE
         ZFLX_SDUST(JL,JAER)=ZFLX_SDUST(JL,JAER)*ZDUEMPOT(JL,JAER)
       ENDIF
       IF (LAERDUSTSIZEVAR) THEN
         IF (JAER <=2 .AND. PDSZ(JL) > 0._JPRB) THEN
           ZFLX_SDUST(JL,JAER)=ZFLX_SDUST(JL,JAER)*PDSZ(JL)
         ELSEIF (PDSZ(JL) > 0._JPRB) THEN
           ZFLX_SDUST(JL,JAER)=ZFLX_SDUST(JL,JAER)-(PDSZ(JL) - 1._JPRB)*(ZFLX_SDUST(JL,1)+ZFLX_SDUST(JL,2))
         ENDIF
       ENDIF

       PAERFLX(JL,JAER)=ZFLX_SDUST(JL,JAER)

     ENDIF
    ENDDO
  ENDDO ! enddo NTYPEAER(2)
 ENDIF ! case NDDUST == 2 or 4

 IF (NDDUST == 3 ) THEN ! case ECMWF formulation

 !- ECMWF dust emission fluxes come in either 3- or 10-size bins
 ! 0.03 - 0.55 - 0.9 - 20.
 ! 0.03 - 0.06 - 0.12 - 0.24 - 0.48 - 0.96 - 1.92 - 3.84 - 7.68 - 15.36 - 30.72

 !!-- for potential dust sources, select land points, snow-free, and zero ice, no wet skin cover
 !!   with fraction of bare soil > 10%, no high vegetation, possible low vegetation < 50% but 
 !!   with soil moisture below moisture corresponding to twice the wilting point (0.171), and 
 !!   a flatish surface (st.dev.orog < 50) with total albedo < 50%

 !-- for potential dust sources, select land points, snow-free, and zero ice, 
 !   no wet skin cover, with some fraction of bare soil, with test on soil 
 !   moisture, and a flatish surface (st.dev.orog < 50) with total albedo < 52%

!DEC$ IVDEP
  DO JL=KIDIA,KFDIA    
!-- default values for non-land points
    LLPDUSTS(JL)=.FALSE.
    PAERMAP(JL,5)=0.0_JPRB
    ZSCC2(JL)=0._JPRB
    ZDEP2(JL)=0._JPRB
    ZLTS2(JL)=0._JPRB  
    IF (PLSM(JL) >= 0.99_JPRB) THEN
      ZREFSPD = PAERMAP(JL,3)
      ZREFRAD = PAERMAP(JL,4)
      ZRADREF = ZREFSPD * ZREFRAD**0.25_JPRB
!-- default min and max of LTS correspond to PWS1 = ZRWPWP and PSW1 = ZRWSAT
      ZLTSMIN(JL) = 0.6_JPRB * ZRADREF       ! ZFSWET = 0.6 
      ZLTSMAX(JL) = 1.2_JPRB * ZRADREF       ! ZFSWET = 1.2
      ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB, (PWS1(JL)-ZRWPWP)/(ZRWSAT-ZRWPWP) ) )
      ZFSWET = 1.2_JPRB+0.2_JPRB*LOG10(ZSWETN)
!-- background lifting threshold speed (defined for all land points)
      ZLTS2(JL) = MIN( ZLTSMAX(JL), MAX( ZLTSMIN(JL), ZFSWET * ZRADREF ))
!-- replace  by simpler test on:
!     absence of snow
!     flatish surface
!     total albedo < 0.52 (no permanent ice)
!     type 8 fraction bare soil > 10%
!     type 4 cover by dry snow-free low vegetated < 50%
!     all other cover types < 0.1%
      IF (PSNS(JL) < ZEPSSNO .AND. PHSDFOR(JL) <= 50._JPRB .AND. PALB(JL) < 0.52_JPRB .AND.&
        & PFRTI(JL,2) < ZEPSARE .AND. PFRTI(JL,3) < ZEPSARE .AND.&! no ice, no wet skin
        & PFRTI(JL,5) < ZEPSARE .AND. PFRTI(JL,6) < ZEPSARE .AND.&! no snow under bare soil-low veg, no dry high veg
        & PFRTI(JL,7) < ZEPSARE .AND.&! no snow under high veg
        & PFRTI(JL,8) > 0.1_JPRB .AND. PFRTI(JL,4) < 0.5_JPRB ) THEN
        LLPDUSTS(JL)=.TRUE.
        PAERMAP(JL,5)=RAERDUB * ZDUEMPOT(JL,1)                         ! for diagnostics only
      ENDIF
    ENDIF
  ENDDO   

! ZFLX_SDUST is positive in kg m-2 s-1
! but ECMWF conventions have PCFLX as a negative upward flux
! input parameters from climatology are:
!  -- soil clay content        (%)
!     dust emission potential  (kg s2 m-5)
!     lifting thereshold speed (m s-1)
!     
  DO JAER=1,NTYPAER(2)

!-- surface source of dust is assumed if LLPTDUSTS is true, and/or UVis albedo > 0.11
!                                        10m wind > threshold = f(soil wetness, mean particle radius)
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
      LLDUST(JL)=.FALSE.
      ZFLX_SDUST(JL,JAER)=0._JPRB

!---------------------------------------------------------------------
!-- ECMWF formulation

      IF (LLPDUSTS(JL)) THEN
        ZDEP2(JL)= RAERDUB * ZDUEMPOT(JL,JAER)
        ZSCC2(JL)= 20._JPRB

!-- Present formulation in MACC (June'11, still kept June'13)
!      use a formula of threshold wind velocity modified from Ginoux et al., 2001
!      based on 1st layer soil wetness and an averaged particle radius
!--    All of the above limes computed above
        PAERLIF(JL,JAER)=ZLTS2(JL)    ! for diagnostics only

!- compare surface 10-m wind with threshold wind velocity

        ZWND3(JL) = MAX(0._JPRB, (ZWNDDU(JL)-ZLTS2(JL)) *ZWNDDU(JL)*ZWNDDU(JL) )

!- preferred approach: flux is based on MODIS-derived UVis_Alb (0.3-0.7 um)
        IF (NALBEDOSCHEME > 0) THEN
!          IF (PALBD(JL,1) >= 0.20_JPRB .AND. PALBD(JL,1) < 0.55_JPRB ) THEN
          IF (PALBD(JL,1) >= 0.08_JPRB .AND. PALBD(JL,1) < 0.55_JPRB ) THEN
            ZFLX_SDUST(JL,JAER)= ZDEP2(JL) * PALBD(JL,1) * ZWND3(JL)
          ENDIF
!-- alternate approach, if MODIS-derived albedo not available, use total albedo
        ELSE 
          ZFLX_SDUST(JL,JAER)= ZDEP2(JL) * PALB(JL) * ZWND3(JL)
        ENDIF

        LLDUST(JL)=.TRUE.
        PAERFLX(JL,JAER) = ZFLX_SDUST(JL,JAER)
      ENDIF
    ENDDO
  ENDDO
ENDIF  ! case NDDUST == 3

!-- PCFLX in kg m-2 s-1

 IF (NDDUST /= 0) THEN
  DO JAER=1,NTYPAER(2)
    INBAER=INBAER+1
    DO JL=KIDIA,KFDIA
      PCFLX(JL,KAERO(INBAER))=-ZFLX_SDUST(JL,JAER) * 1.E+00_JPRB
    ENDDO


!-- if no vertical diffusion, distribute the flux in layers with scale height
!-- between half-levels IHTST-1 and KLEV
    IF (.NOT.LVDFTRAC) THEN
      DO JK=IHTST,KLEV
        DO JL=KIDIA,KFDIA
          ZDETAH(JL,JK) = ZETAH(JL,JK)**ZHDD - ZETAH(JL,JK-1)**ZHDD
          PTENC(JL,JK,KAERO(INBAER)) = PTENC(JL,JK,KAERO(INBAER))+PCFLX(JL,KAERO(INBAER))*ZDETAH(JL,JK)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
 ENDIF
ENDIF

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
  CALL SATUR (KIDIA, KFDIA, KLON  , KTDIA, KLEV, YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
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

