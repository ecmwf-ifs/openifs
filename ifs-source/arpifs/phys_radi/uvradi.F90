! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UVRADI &
  &( YDMODEL,KIDIA, KFDIA, KLON, KLEV , KAER , KSW  , KUV  , &
  &  PAERO, PALBD, PALBP, PAPH, PAP  , PCCNL, PCCNO, PCLDF, &
  &  PGELAM,PCLON, PSLON, PDP , PGEMU, PMU0 , POZ  , PQ   , PQSAT , &
  &  PQICE, PQLI , PLSM, PT   , PTL  , PUVC , PUVT , &
  &  PUVCTOT,PUVTTOT, VUVP1)

!**** *UVRADI* - COMPUTES THE ULTRA-VIOLET SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE UV SURFACE RADIATION FLUXES IN N
!          SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980)
!          OR RRTMG-SW

!**   INTERFACE.
!     ----------

!          *UVRADI* IS CALLED FROM *CALLPAR*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
! PAERO  (KLON,KLEV,NACTAERO): PROGNOSTIC AEROSOLS
! PALBD  (KLON,KSW)    : DIFFUSE SURFACE ALBEDO
! PALBP  (KLON,KSW)    : DIRECT/PARALLEL SURFACE ALBEDO
! PAPH   (KLON,KLEV+1) : HALF-LEVEL PRESSURE (Pa)
! PAP    (KLON,KLEV)   : FULL-LEVEL PRESSURE (Pa)
! PCCNL  (KLON)        : NUMBER OF CCN OVER LAND  (m-3)
! PCCNO  (KLON)        : NUMBER OF CCN OVER OCEAN (m-3)
! PCLDF  (KLON,KLEV)   : CLOUD FRACTION
! PGELAM (KLON)        : LONGITUDE (radian)
! PCLON  (KLON)        : COSINE OF LONGITUDE
! PSLON  (KLON)        : SINE OF LONGITUDE
! PDP    (KLON,KLEV)   : PRESSURE THICKNESS OF LAYER (Pa)
! PGEMU  (KLON)        : LATITUDE (radian)
! PMU0   (KLON)        : COSINE SOLAR ZENITH ANGLE
! POZ    (KLON,KLEV)   : OZONE MIXING RATIO (kg/kg) -- prognostic --
! PQ     (KLON,KLEV)   : HUMIDITY (kg/kg)
! PQSAT  (KLON,KLEV)   : HUMIDITY AT SATURATION (kg/kg)
! PQICE  (KLON,KLEV)   : CLOUD ICE MIXING RATIO (kg/kg)
! PQLI   (KLON,KLEV)   : CLOUD LIQUID WATER M.RATIO (kg/kg)
! PLSM   (KLON)        : LAND-SEA MASK
! PT     (KLON,KLEV)   : TEMPERATURE (K)
! PTL    (KLON)        : SURFACE/SKIN TEMPERATURE (K)

!     ==== OUTPUTS ===
! PUVC   (KLON,KUV)    : SURFACE UV CLEAR-SKY SPECTRAL FLUX (W m-2 nm-1)
! PUVT   (KLON,KUV)    : SURFACE UV TOTAL SKY SPECTRAL FLUX (W m-2 nm-1)
! PUVCTOT(KLON)        : SURFACE UV CLEAR-SKY FLUX (W m-2)
! PUVTTOT(KLON)        : SURFACE UV TOTAL SKY FLUX (W m-2)
! PUVPOUT(KLON,KLEV)   : OUTPUTS RADIANCES

!     METHOD.
!     -------
!          trimmed down version of SW scheme, but with increased spectral 
!          resolution in the UV

!          1. COMPUTES FLUXES IN U.V. INTERVALS         (UVFLX)

!     EXTERNALS.
!     ----------

!          *UVFLX*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        ORIGINAL : 2005-10-04 based on *RADLSW* and *SW* but top to bottom

!     MODIFICATIONS.
!     --------------
!     N.Semane+P.Bechtold 04-10-2012 add RPLRG for small planet
!     K. Yessad (July 2014): Move some variables.
!     A.Bozzo (Oct 2014) revision based on RRTMG-SW v3.9
!     A.Bozzo (Jun 2016) fixed bug in the aerosol optical properties for prognostic types
!     A.Bozzo (Sep 2016) fixed bug in cloud overlap assumptions
!     ------------------------------------------------------------------

USE TYPE_MODEL     , ONLY : MODEL
USE PARKIND1       , ONLY : JPIM, JPRB
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST         , ONLY : RG, RD, RMD, RMO3, RNAVO, RTT, RPI
USE YOMCT0         , ONLY : LIFSMIN, LIFSTRAJ
USE YOMCT3         , ONLY : NSTEP
USE YOERDU         , ONLY : REPLOG, REPSC, REPSCT, REPSCW
USE YOETHF         , ONLY : RTICE
USE YOMDYNCORE     , ONLY : RPLRG
USE RRTMG_SW_SPCVRT, ONLY : SPCVRT_SW
USE YOMLUN         , ONLY : NULOUT
USE VARIABLE_MODULE, ONLY : VARIABLE_3D

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA

INTEGER(KIND=JPIM),INTENT(IN)    :: KAER
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KUV

REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCNL(KLON),PCCNO(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON),PCLON(KLON),PSLON(KLON),PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON),PMU0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQICE(KLON,KLEV),PQLI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV),PTL(KLON)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVC(KLON,KUV), PUVT(KLON,KUV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVCTOT(KLON) , PUVTTOT(KLON) 
TYPE(VARIABLE_3D) ,INTENT(INOUT) :: VUVP1

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZAKI(KLON,2) &
  &, ZALBD(KLON)         , ZALBP(KLON)      , ZC1J(KLON,KLEV+1) &
  &, ZCLD(KLON,KLEV)     , ZCLEAR(KLON)     , ZCLOUD(KLON)      &
  &, ZDSIG(KLON,KLEV)    , ZFACT(KLON)      &
  &, ZCD(KLON,KUV,KLEV+1), ZFD(KLON,KUV,KLEV+1) &
  &, ZCDUV(KLON,KLEV+1)  , ZFDUV(KLON,KLEV+1)   & 
  &, ZQAER(KLON,6,KLEV)  , ZQOZ(KLON,KLEV)  , ZQOZC(KLON,KLEV)  &
  &, ZRMU(KLON),ZSEC(KLON) &
  &, ZSIGO(KLON)         , ZSIGN(KLON)      , ZTH(KLON,KLEV+1)  &
  &, ZTAU(KLON,KLEV)     , ZCG(KLON,KLEV)   , ZOMEGA(KLON,KLEV) &
  &, ZFIWP(KLON,KLEV)    , ZFLWP(KLON,KLEV) &
  &, ZIWC(KLON)          , ZLWC(KLON)       &
  &, ZCGAK(KLON,KLEV)    , ZPIZAK(KLON,KLEV), ZTAUAK(KLON,KLEV) &
  &, ZCGAZ(KLON,KLEV)    , ZPIZAZ(KLON,KLEV), ZTAUAZ(KLON,KLEV) &
  &, ZDESR(KLON,KLEV)    , ZRADIP(KLON,KLEV), ZRADLP(KLON,KLEV) &
  &, ZCOZ(KLON)          , ZRAYL(KLON)      &
  &, ZUD(KLON,2,KLEV+1)  , ZUDAC(KLON,2)    &
  &, ZGAER(KLON,8,KLEV)  , ZOAER(KLON,8,KLEV)&
  &, ZDUM(KLON),ZTAUO3(KLON,KLEV) &
  &, ZTRAY(KLON,KLEV), ZTAUAER(KLON,8,KLEV) &
  &, ZDUMAER(KLON,KLEV,12)

REAL(KIND=JPRB) :: ZTCRAY(KLON,KUV),ZTCO3(KLON,KUV),ZTCAER(KLON,KUV)
REAL(KIND=JPRB) :: ZUVTSI(KLON),ZUVDDCTOT(KLON),ZUVDDTTOT(KLON)
REAL(KIND=JPRB) :: ZUVDDC(KLON,KUV), ZUVDDT(KLON,KUV) 

REAL(KIND=JPRB) :: Z1S, Z2S, Z3S, ZCRAE

REAL(KIND=JPRB) :: ZASYMX, ZGI   , ZGL   , ZIWGKG, ZLWGKG,&
 & ZOI   , ZOL   , ZOMGMX, ZTAUMX, ZTEMPC, Z1RADI,&
 & ZTOI, ZTOL, ZDPOG, ZPODT  

REAL(KIND=JPRB) :: ZALND  , ZASEA  , ZD     , ZDEN  , ZNTOT, ZNUM  , ZRATIO, &
 & ZBETAI, ZOMGI , ZTCELS , ZFSR , ZAIWC , &
 & ZBIWC , &
 & ZDEFRE, ZREFDE, ZALPHA1, &
 & ZCSTO , ZCSTR , ZDENB  , ZEFACT , ZFF    , ZGAR  , ZT250, ZTFACT

!output arrays from rrtm-sw
REAL(KIND=JPRB) :: ZBBFD(KLON,KLEV+1,KUV), ZBBFU(KLON,KLEV+1,KUV), &
 & ZBBCD(KLON,KLEV+1,KUV), ZBBCU(KLON,KLEV+1,KUV), &
 & ZBBFDDIR(KLON,KLEV+1,KUV), ZBBCDDIR(KLON,KLEV+1,KUV) 

INTEGER(KIND=JPIM) :: INWAVL, ITWAVL(19),IOVLP_UV
INTEGER(KIND=JPIM) :: JAER, JAEUV, JALUV, JCLUV, JK, JL, JPUV, JUV

REAL(KIND=JPRB) :: ZUVPOUT(KLON,KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "aer_rad.intfb.h"
#include "radaca.intfb.h"
#include "uvflxa.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UVRADI',0,ZHOOK_HANDLE)
ASSOCIATE(YDEUVRAD=>YDMODEL%YRML_PHY_RAD%YREUVRAD,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDEOVLP=>YDMODEL%YRML_PHY_RAD%YREOVLP,YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, &
 & YDECLD=>YDMODEL%YRML_PHY_EC%YRECLD,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI,YDEAERD=>YDMODEL%YRML_PHY_RAD%YREAERD, &
 & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
 & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NDIM=>YGFL%NDIM, YUVP=>YGFL%YUVP, &
 & LUVPOUT=>YGFL%LUVPOUT, &
 & RCAEROS=>YDEAERD%RCAEROS, &
 & REPSEC=>YDECLD%REPSEC, &
 & RA1OVLP=>YDEOVLP%RA1OVLP, &
 & LEO3CH=>YDEPHY%LEO3CH, &
 & KMODTS=>YDERAD%KMODTS, LCCNL=>YDERAD%LCCNL, LCCNO=>YDERAD%LCCNO, &
 & NOVLP=>YDERAD%NOVLP, NRADIP=>YDERAD%NRADIP, NRADLP=>YDERAD%NRADLP, &
 & RCCNLND=>YDERAD%RCCNLND, RCCNSEA=>YDERAD%RCCNSEA, &
 & RRAE=>YDERDI%RRAE, &
 & IPUV=>YDEUVRAD%IPUV, JCOP=>YDEUVRAD%JCOP, LUVAERP=>YDEUVRAD%LUVAERP, &
 & LUVTDEP=>YDEUVRAD%LUVTDEP, RASA=>YDEUVRAD%RASA, RASB=>YDEUVRAD%RASB, &
 & RASC=>YDEUVRAD%RASC, RASD=>YDEUVRAD%RASD, RASE=>YDEUVRAD%RASE, &
 & RASF=>YDEUVRAD%RASF, RAYUVB=>YDEUVRAD%RAYUVB, RCGUVA=>YDEUVRAD%RCGUVA, &
 & RCIEAS=>YDEUVRAD%RCIEAS, RFA0=>YDEUVRAD%RFA0, RFA1=>YDEUVRAD%RFA1, &
 & RFB0=>YDEUVRAD%RFB0, RFB1=>YDEUVRAD%RFB1, RFB2=>YDEUVRAD%RFB2, &
 & RFB3=>YDEUVRAD%RFB3, RFC0=>YDEUVRAD%RFC0, RFC1=>YDEUVRAD%RFC1, &
 & RFC2=>YDEUVRAD%RFC2, RFC3=>YDEUVRAD%RFC3, RFCAER=>YDEUVRAD%RFCAER, &
 & RK250=>YDEUVRAD%RK250, RMUZUV=>YDEUVRAD%RMUZUV, RPIUVA=>YDEUVRAD%RPIUVA, &
 & RSUVB=>YDEUVRAD%RSUVB, RTAUVA=>YDEUVRAD%RTAUVA, RTUV1=>YDEUVRAD%RTUV1, &
 & RTUV2=>YDEUVRAD%RTUV2, RUVLAM=>YDEUVRAD%RUVLAM, &
 & NSTART=>YDRIP%NSTART)
!     ------------------------------------------------------------------

ZUVPOUT(:,:)=0.0_JPRB

!*         1.     USEFUL CLOUD-RELATED QUANTITIES
!                 -------------------------------

!cloud overlap assumption in the processor
!0=only cloud fraction used
!1=max-random 
!2=maximum
!3=random
!4=exp-max-rand (coeffs defined only up to 60 levels?)
IOVLP_UV=1  
!WRITE(UNIT=NULOUT,FMT='('' IOVLP_UV = '',I1)') IOVLP_UV


ZCRAE=RRAE*(RRAE+2.0_JPRB)

DO JL=KIDIA,KFDIA
  IF (PMU0(JL) > 1.E-10_JPRB) THEN
    ZRMU(JL)=RRAE/(SQRT(PMU0(JL) * PMU0(JL) + ZCRAE)-PMU0(JL))
  ELSE
    ZRMU(JL)=RRAE/SQRT(ZCRAE)
  ENDIF
  ZSEC(JL)=1.0_JPRB/ZRMU(JL)
  ZC1J(JL,1)=0.0_JPRB
  ZSIGO(JL) = PAPH(JL,1)
  ZCLEAR(JL)=1.0_JPRB
  ZCLOUD(JL)=0.0_JPRB
  ZUDAC(JL,1)=0.0_JPRB
  ZUDAC(JL,2)=0.0_JPRB
ENDDO
DO JK = 1 , KLEV
  ZALPHA1=RA1OVLP(JK)  
  DO JL = KIDIA,KFDIA
    ZSIGN(JL) = PAPH(JL,JK+1)
    ZDSIG(JL,JK) = (ZSIGN(JL) - ZSIGO(JL))/PAPH(JL,KLEV+1)
    ZSIGO(JL) = ZSIGN(JL)

    IF (IOVLP_UV == 1) THEN !max-ran
      ZCLEAR(JL)=ZCLEAR(JL)&
       & *(1.0_JPRB-MAX(PCLDF(JL,JK),ZCLOUD(JL)))&
       & /(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC))  
      ZC1J(JL,JK+1)= 1.0_JPRB - ZCLEAR(JL)
      ZCLOUD(JL) = PCLDF(JL,JK)
    ELSEIF (IOVLP_UV == 2) THEN !maximum
      ZCLOUD(JL) = MAX(PCLDF(JL,JK),ZCLOUD(JL))
      ZC1J(JL,JK+1) = ZCLOUD(JL)
    ELSEIF (IOVLP_UV == 3) THEN !random
      ZCLEAR(JL) = ZCLEAR(JL)*(1.0_JPRB-PCLDF(JL,JK))
      ZCLOUD(JL) = 1.0_JPRB - ZCLEAR(JL)
      ZC1J(JL,JK+1) = ZCLOUD(JL)
    ELSEIF (IOVLP_UV == 4) THEN !exp-max-rand
!** Hogan & Illingworth (2001)      
      ZCLEAR(JL)=ZCLEAR(JL)*( &
       & ZALPHA1*(1.0_JPRB-MAX(PCLDF(JL,JK),ZCLOUD(JL))) &
       & /(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC)) &
       & +(1.0_JPRB-ZALPHA1)*(1.0_JPRB-PCLDF(JL,JK)) )  
      ZC1J(JL,JK+1) = 1.0_JPRB - ZCLEAR(JL) 
      ZCLOUD(JL) = PCLDF(JL,JK)
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZCLEAR(JL)=1.0_JPRB-ZC1J(JL,KLEV+1)
!-- no fiddling -------------------------
  IF (IOVLP_UV == 0) ZCLEAR(JL)=0._JPRB
ENDDO
JL=KIDIA
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (ZCLEAR(JL) < 1.0_JPRB) THEN
      ZCLD(JL,JK)=PCLDF(JL,JK)/(1.0_JPRB-ZCLEAR(JL))
    ELSE
      ZCLD(JL,JK)=0.0_JPRB
    ENDIF
    ZCLD(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZCLD(JL,JK)))
!-- no fiddling -------------------------------------
!    ZCLD(JL,JK)=PCLDF(JL,JK)
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         2.     OZONE AND AEROSOLS
!                 ------------------

!-- BY DEFAULT, USE TEGEN ET AL. CLIMATOLOGY OF AEROSOLS

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    Z1S = PAP(JL,JK-1)*(PAP(JL,JK)-PAPH(JL,JK))
    Z2S = PAP(JL,JK)  *(PAPH(JL,JK)-PAP(JL,JK-1))
    Z3S = PAPH(JL,JK) *(PAP(JL,JK)-PAP(JL,JK-1))
    ZTH(JL,JK)=(PT(JL,JK-1)*Z1S+PT(JL,JK)*Z2S)/Z3S
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZTH(JL,1) =&
   & PT(JL,1)-PAP(JL,1) *(PT(JL,1)-ZTH(JL,2)) / (PAP(JL,1)-PAPH(JL,2))    
  ZTH(JL,KLEV+1)= PTL(JL)
ENDDO

!-- diagnostic aerosols from climatologies by default

CALL  RADACA ( YDEAERD,YDERAD, YDRIP, KIDIA  , KFDIA, KLON , KLEV,&
  & PAPH , PGELAM, PGEMU, PCLON, PSLON, ZTH , &
  & ZQAER, ZDUMAER, ZQOZC, ZDUM,  &
  & ZDUM, ZDUM, ZDUM, ZDUM, ZDUM, ZDUM ) 

DO JAER = 1 , 6
  DO JK = 1 , KLEV
    DO JL = KIDIA,KFDIA
      ZQAER(JL,JAER,JK) = ZQAER(JL,JAER,JK)*RFCAER
    ENDDO
  ENDDO
ENDDO

!-- prognostic aerosols formatted as the diagnostic aerosols

IF (NACTAERO > 0 .AND. LUVAERP) THEN

  INWAVL=1

  ITWAVL(:)=0
  ITWAVL(1)=9

! -- profile of aerosol optical parameters at 550 nm

  CALL AER_RAD(YDEAERATM, YDMODEL%YRML_PHY_AER, YDCOMPO, YGFL, &
   & KIDIA, KFDIA, KLON, 1, KLEV , NACTAERO, INWAVL, ITWAVL, NSTART, NSTEP, &
   & PAP, PAPH, PQ, PT, PAERO, &
   & ZTAUAER, ZOAER, ZGAER )

  ! NB: ZTAUAER have SS, DD, OM, BC (inc volcanic ash), SU (inc volcanic), NI, AM in that order and have to be
  !     redistributed over the old Tegen "radiative" aerosols
  !     1 = land   = organic + sulphate + nitrate + ammonium
  !     2 = sea    = sea salt
  !     3 = desert = dust
  !     4 = urban  = black carbon
  !     5 = volcanic
  !     6 = stratospheric background
  !  In present configuration of prognostic aerosols, types 5 (volcanic =0.)
  !     and 6 (stratospheric background) are still given by climatology.

  DO JK=1,KLEV
     DO JL=KIDIA,KFDIA
        ZQAER(JL,1,JK)=ZTAUAER(JL,3,JK)+ZTAUAER(JL,5,JK)+ZTAUAER(JL,6,JK)+ZTAUAER(JL,7,JK) ! organic+sulphate+nitrate+ammonium
        ZQAER(JL,2,JK)=ZTAUAER(JL,1,JK)                  !  sea salt
        ZQAER(JL,3,JK)=ZTAUAER(JL,2,JK)                   ! desert dust
        ZQAER(JL,4,JK)=ZTAUAER(JL,4,JK)                   ! black carbon
     ENDDO
  ENDDO

ENDIF

!-- NB: the UV processor uses the prognostic ozone POZ if it is active

ZCSTO = 1.E-20_JPRB*RNAVO / (RMO3*(RG/RPLRG)*10._JPRB)

! factor 1.E-20 as O3 cross-section does not include this factor
DO JL=KIDIA,KFDIA
  ZCOZ(JL)=0._JPRB
ENDDO
!-- LEO3CH = .T. UV diagnostics from prognostic O3, .F. from climatological O3
!depending from the choice in uvradi_layer this can come from the full chemistry
!or the parametrized one
IF (LEO3CH) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZQOZ(JL,JK)=POZ(JL,JK)*PDP(JL,JK)*ZCSTO
      ZCOZ(JL)   =ZCOZ(JL) + ZQOZ(JL,JK)
    ENDDO
  ENDDO
ELSE
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZQOZ(JL,JK)=ZQOZC(JL,JK)*ZCSTO
      ZCOZ(JL)   =ZCOZ(JL) + ZQOZ(JL,JK)
    ENDDO
  ENDDO
ENDIF

!print *,'KAER=',KAER

IF (KAER == 0) THEN
  DO JAER = 1 , 6
    DO JK = 1 , KLEV
      DO JL = KIDIA,KFDIA
        ZQAER(JL,JAER,JK) = RCAEROS
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JK = 1 , KLEV
    DO JL = KIDIA,KFDIA
      ZQAER(JL,5,JK) = RCAEROS
    ENDDO
  ENDDO
ENDIF



!     ------------------------------------------------------------------

!*         3.     RAYLEIGH, CLOUD AND AEROSOLS OPTICAL PROPERTIES 
!                 -----------------------------------------------

ZREFDE = 0.64952_JPRB
ZDEFRE = 1.0_JPRB / ZREFDE

DO JUV=1,KUV
  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFD(JL,JUV,JK) =0.0_JPRB
      ZCD(JL,JUV,JK) =0.0_JPRB
      ZBBFD(JL,JK,JUV) = 0.0_JPRB
      ZBBFU(JL,JK,JUV) = 0.0_JPRB
      ZBBCD(JL,JK,JUV) = 0.0_JPRB
      ZBBCU(JL,JK,JUV) = 0.0_JPRB
      ZBBFDDIR(JL,JK,JUV) = 0.0_JPRB
      ZBBCDDIR(JL,JK,JUV) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

!*         3.1    PREPARATORY FOR CLOUD OPTICAL PROPERTIES
!                 ----------------------------------------

! OPTICAL PROPERTIES FOR LIQUID WATER CLOUDS FROM SLINGO (1989) 
! AND FOR ICE WATER CLOUDS FROM FU (1996)

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

! --- LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
    IF (PCLDF(JL,JK) > REPSC ) THEN
      ZLWGKG=MAX(PQLI(JL,JK)*1000.0_JPRB,0.0_JPRB)
      ZIWGKG=MAX(PQICE(JL,JK)*1000.0_JPRB,0.0_JPRB)
      ZLWGKG=ZLWGKG/PCLDF(JL,JK)
      ZIWGKG=ZIWGKG/PCLDF(JL,JK)
    ELSE
      ZLWGKG=0.0_JPRB
      ZIWGKG=0.0_JPRB
    ENDIF
    ZDPOG=PDP(JL,JK)/(RG/RPLRG)
    ZFLWP(JL,JK)= ZLWGKG*ZDPOG
    ZFIWP(JL,JK)= ZIWGKG*ZDPOG
    ZPODT=PAP(JL,JK)/(RD*PT(JL,JK))
    ZLWC(JL)=ZLWGKG*ZPODT
    ZIWC(JL)=ZIWGKG*ZPODT
  ENDDO

!*         3.2    CLOUD OPTICAL PROPERTIES: EFFECTIVE DROPLET RADIUS AND PARTICLE SIZE
!                 --------------------------------------------------------------------

  DO JL = KIDIA,KFDIA
!---- EFFECTIVE RADIUS FOR WATER, ICE AND RAIN PARTICLES

! very old parametrization as f(pressure)

    IF (NRADLP == 0) THEN
!-- very old parametrization as f(pressure) ERA-15
      ZRADLP(JL,JK)=10.0_JPRB + (100000.0_JPRB-PAP(JL,JK))*3.5_JPRB

    ELSEIF (NRADLP == 1) THEN
! simple distinction between land (10) and ocean (13) Zhang and Rossow
      IF (PLSM(JL) < 0.5_JPRB) THEN
        ZRADLP(JL,JK)=13.0_JPRB
      ELSE
        ZRADLP(JL,JK)=10.0_JPRB
      ENDIF
      
    ELSEIF (NRADLP >= 2) THEN
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
        ZRADLP(JL,JK)=100.0_JPRB*EXP(0.333_JPRB*LOG(ZNUM/ZDEN))
        ZRADLP(JL,JK)=MAX(ZRADLP(JL,JK), 4.0_JPRB)
        ZRADLP(JL,JK)=MIN(ZRADLP(JL,JK),16.0_JPRB)
      ELSE
        ZRADLP(JL,JK)=4.0_JPRB
      ENDIF
    ENDIF  
  ENDDO

  DO JL = KIDIA,KFDIA

! diagnosing the ice particle effective radius/diameter

!- ice particle effective radius =f(T) from Liou and Ou (1994)
 
    IF (PT(JL,JK) < RTICE) THEN
      ZTEMPC=PT(JL,JK)-RTT
    ELSE
      ZTEMPC=RTICE-RTT
    ENDIF
    ZRADIP(JL,JK)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
      & 0.0012_JPRB))    

    IF (NRADIP == 0) THEN
!-- fixed 40 micron effective radius
      ZRADIP(JL,JK)= 40.0_JPRB
      ZDESR(JL,JK) = ZDEFRE * ZRADIP(JL,JK)
      
    ELSEIF (NRADIP == 1) THEN 

!-- old formulation based on Liou & Ou (1994) temperature (40-130microns)    
      ZRADIP(JL,JK)=MAX(ZRADIP(JL,JK),40.0_JPRB)
      ZDESR(JL,JK) = ZDEFRE * ZRADIP(JL,JK)
      
    ELSEIF (NRADIP == 2) THEN  
!-- formulation following Jakob, Klein modifications to ice content    
      ZRADIP(JL,JK)=MAX(ZRADIP(JL,JK),30.0_JPRB)
      ZRADIP(JL,JK)=MIN(ZRADIP(JL,JK),60.0_JPRB)
      ZDESR(JL,JK)= ZDEFRE * ZRADIP(JL,JK)
 
    ELSEIF (NRADIP == 3  ) THEN
 
!- ice particle effective radius =f(T,IWC) from Sun and Rikus (1999)
! revised by Sun (2001)
      IF (ZIWC(JL) > 0.0_JPRB ) THEN
        ZTEMPC = PT(JL,JK)-83.15_JPRB
        ZTCELS = PT(JL,JK)-RTT
        ZFSR = 1.2351_JPRB +0.0105_JPRB * ZTCELS
! Sun, 2001 (corrected from Sun & Rikus, 1999)
        ZAIWC = 45.8966_JPRB * ZIWC(JL)**0.2214_JPRB
        ZBIWC = 0.7957_JPRB * ZIWC(JL)**0.2535_JPRB
        ZDESR(JL,JK) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)
        ZDESR(JL,JK) = MIN ( MAX( ZDESR(JL,JK), 30.0_JPRB), 155.0_JPRB)
        ZRADIP(JL,JK)= ZREFDE * ZDESR(JL,JK)
      ELSE
        ZDESR(JL,JK) = 92.5_JPRB
        ZRADIP(JL,JK)= ZREFDE * ZDESR(JL,JK)
      ENDIF  
    ENDIF  
    
  ENDDO
ENDDO


!*         3.3    UV-SW OPTICAL PROPERTIES: RAYLEIGH, CLOUDS AND AEROSOLS
!                 -------------------------------------------------------

DO JPUV=1,KUV
  JUV  =IPUV(JPUV)
  JCLUV=JCOP(JUV) 
  JAEUV=JCOP(JUV)
  JALUV=1

!*         3.3.1  OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------

! factor 10 to convert Rayleigh scattering cross section (cm2) to optical thickness (N.D.)
  ZCSTR = RAYUVB(JUV) * RNAVO / (RMD*(RG/RPLRG)*10._JPRB)

  DO JL = KIDIA,KFDIA
    IF (PMU0(JL) >= RMUZUV) THEN 
      ZFACT(JL) = PMU0(JL) * RSUVB(JUV)
    ELSE
      ZFACT(JL) = 0._JPRB
    ENDIF
    ZRAYL(JL) = ZCSTR * PAPH(JL,KLEV+1) 
    ZTCRAY(JL,JPUV)=ZRAYL(JL)
    ZTCAER(JL,JPUV)=0.0_JPRB
  ENDDO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZALBD(JL)=PALBD(JL,JALUV)
      ZALBP(JL)=PALBP(JL,JALUV)
      ZTAU(JL,JK)  = 0.0_JPRB
      ZOMEGA(JL,JK)= 1.0_JPRB
      ZCG(JL,JK)   = 0.0_JPRB
      ZCGAZ(JL,JK) = 0.0_JPRB
      ZPIZAZ(JL,JK)= 0.0_JPRB
      ZTAUAZ(JL,JK)= 0.0_JPRB
    ENDDO


!*         3.3.2  OPTICAL THICKNESS FOR MIX OF RAYLEIGH AND AEROSOLS LAYER-BY-LAYER
!                 -----------------------------------------------------------------

    DO JAER=1,6
      DO JL = KIDIA,KFDIA
        ZTAUAZ(JL,JK)=ZTAUAZ(JL,JK)+ZQAER(JL,JAER,JK)*RTAUVA(JAEUV,JAER)
        ZPIZAZ(JL,JK)=ZPIZAZ(JL,JK)+ZQAER(JL,JAER,JK)&
         & * RTAUVA(JAEUV,JAER)*RPIUVA(JAEUV,JAER)  
        ZCGAZ(JL,JK) = ZCGAZ(JL,JK)+ZQAER(JL,JAER,JK)&
         & * RTAUVA(JAEUV,JAER)*RPIUVA(JAEUV,JAER)*RCGUVA(JAEUV,JAER)  
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      ZTCAER(JL,JPUV)=ZTCAER(JL,JPUV)+ZTAUAZ(JL,JK)
    ENDDO

    DO JL = KIDIA,KFDIA
      IF (KAER /= 0) THEN
        ZCGAK(JL,JK) = ZCGAZ(JL,JK)/ZPIZAZ(JL,JK)
        ZPIZAK(JL,JK)=ZPIZAZ(JL,JK)/ZTAUAZ(JL,JK)
        ZTAUAK(JL,JK)=ZTAUAZ(JL,JK)
        ZTRAY(JL,JK)  = ZRAYL(JL) * ZDSIG(JL,JK)
!!!! wrong  ZRATIO = ZTRAY / (ZTRAY + ZTAUAK(JL,JK))
        ZGAR = ZCGAK(JL,JK)
        ZFF  = ZGAR * ZGAR

        ZDENB = ZTRAY(JL,JK)+ZTAUAK(JL,JK)*(1.0_JPRB-ZPIZAK(JL,JK)*ZFF)
        ZRATIO= ZTRAY(JL,JK) / ZDENB

        ZTAUAZ(JL,JK)= ZTRAY(JL,JK)+ ZTAUAK(JL,JK)*(1.0_JPRB-ZPIZAK(JL,JK)*ZFF)
        ZCGAZ(JL,JK) = ZGAR * (1.0_JPRB - ZRATIO) / (1.0_JPRB + ZGAR)
        ZPIZAZ(JL,JK)= ZRATIO+(1.0_JPRB-ZRATIO)*ZPIZAK(JL,JK)*(1.0_JPRB-ZFF)&
         & / (1.0_JPRB - ZPIZAK(JL,JK) * ZFF)  
      ELSE
        ZTRAY(JL,JK) = ZRAYL(JL) * ZDSIG(JL,JK)
        ZTAUAZ(JL,JK) = ZTRAY(JL,JK)
        ZCGAZ(JL,JK)  = 0.0_JPRB
        ZPIZAZ(JL,JK) = 1.0_JPRB-REPSCT
        ZTAUAK(JL,JK) = 0.0_JPRB
        ZCGAK(JL,JK)  = 0.0_JPRB
        ZPIZAK(JL,JK) = 1.0_JPRB-REPSCT
      ENDIF
    ENDDO


!*         3.3.3  OPTICAL THICKNESS FOR CLOUDS LAYER-BY-LAYER
!                 -------------------------------------------

    DO JL = KIDIA,KFDIA
      ZTOL=0.0_JPRB
      ZGL =0.0_JPRB
      ZOL =0.0_JPRB
      ZTOI=0.0_JPRB
      ZGI =0.0_JPRB
      ZOI =0.0_JPRB
      IF (PCLDF(JL,JK) > REPSC .AND. ZFLWP(JL,JK)+ZFIWP(JL,JK) > REPSCW ) THEN
        IF (ZFLWP(JL,JK) > REPSCW ) THEN
!-- SW: Slingo, 1989
          ZTOL = ZFLWP(JL,JK)*(RASA(JCLUV)+RASB(JCLUV)/ZRADLP(JL,JK))
          ZGL  = RASE(JCLUV)+RASF(JCLUV)*ZRADLP(JL,JK)
          ZOL  = 1._JPRB - RASC(JCLUV)-RASD(JCLUV)*ZRADLP(JL,JK)
        ENDIF

        IF (ZFIWP(JL,JK) > REPSCW ) THEN
!-- SW: Fu 1996
          Z1RADI = 1.0_JPRB / ZDESR(JL,JK)
          ZBETAI = RFA0(JCLUV)+Z1RADI* RFA1(JCLUV)
          ZTOI = ZFIWP(JL,JK) * ZBETAI
          ZOMGI= RFB0(JCLUV)+ZDESR(JL,JK)*(RFB1(JCLUV) + ZDESR(JL,JK) &
           &   *(RFB2(JCLUV)+ZDESR(JL,JK)* RFB3(JCLUV) ))            
          ZOI  = 1.0_JPRB - ZOMGI
          ZGI  = RFC0(JCLUV)+ZDESR(JL,JK)*(RFC1(JCLUV) + ZDESR(JL,JK) &
           &   *(RFC2(JCLUV)+ZDESR(JL,JK)* RFC3(JCLUV) )) 
          ZGI  = MIN(1.0_JPRB, ZGI)
        ENDIF

!  - MIX of WATER and ICE CLOUDS
        ZTAUMX= ZTOL + ZTOI
        ZOMGMX= ZTOL*ZOL + ZTOI*ZOI
        ZASYMX= ZTOL*ZOL*ZGL + ZTOI*ZOI*ZGI

        ZASYMX= ZASYMX/ZOMGMX
        ZOMGMX= ZOMGMX/ZTAUMX

! --- SW FINAL CLOUD OPTICAL PARAMETERS

        ZTAU(JL,JK)  = ZTAUMX
        ZOMEGA(JL,JK)= ZOMGMX
        ZCG(JL,JK)   = ZASYMX
      ENDIF
    ENDDO
  ENDDO


!*         3.3.4  EFFECTIVE ABSORPTION COEFFICIENT FOR O3 AND NO2 INCLUDING TEMPERATURE EFFECT
!                 ----------------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    ZUDAC(JL,1)=0.0_JPRB
    ZUDAC(JL,2)=0.0_JPRB
    ZUD(JL,1,1)=0.0_JPRB
    ZUD(JL,2,1)=0.0_JPRB
  ENDDO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF (LUVTDEP) THEN
        ZT250= PT(JL,JK)-250._JPRB
        ZTFACT=ZT250*( RTUV1(JUV) + (ZT250* RTUV2(JUV) ))
        ZEFACT=EXP(ZTFACT)
      ELSE
        ZEFACT=1.0_JPRB
      ENDIF
      ZUD(JL,1,JK+1)=ZQOZ(JL,JK)*ZEFACT
      ZUDAC(JL,1)=ZUDAC(JL,1)+ZQOZ(JL,JK)*ZSEC(JL)
      ZTAUO3(JL,JK)=RK250(JUV)*ZEFACT*ZQOZ(JL,JK) !O3 optical depth at the wavelength JUV
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
!    ZTR=EXP(-RK250(JUV)*ZUDAC(JL,1))
!    ZAKI(JL,1)=-LOG( ZTR(JL,1) ) / ZUDAC(JL,1)
    ZAKI(JL,1)=MAX(RK250(JUV), 1.E-18_JPRB)
    ZAKI(JL,2)=1.E-12_JPRB
  ENDDO

!     ------------------------------------------------------------------

!*         4.     CALL FLUX ROUTINE 
!                 -----------------

!  CALL UVFLX &
!   & ( KIDIA, KFDIA , KLON  , KLEV , JUV   , &
!   &   ZALBD, ZALBP , ZCG   , ZCLD , ZCLEAR, &
!   &   ZDSIG, ZOMEGA, ZRAYL , ZRMU , ZSEC  , &
!   &   ZTAU , ZCGAZ , ZPIZAZ,ZTAUAZ, ZUD   , &
!   &   ZCDUV, ZFDUV &
!   & ) 


IF (KMODTS == 0) THEN
!Morcrette, Fouquart&Bonnel SW code adapted to UV 

  CALL UVFLXA &
   &( YDMODEL%YRML_PHY_RAD,YDECLD, &
   &  KIDIA, KFDIA , KLON  , KLEV  , JUV  , &
   &  ZAKI , ZALBD , ZALBP , ZCG   , ZCLD , ZCLEAR, &
   &  ZDSIG, ZOMEGA, ZRAYL , ZRMU  , ZSEC , &
   &  ZTAU , ZCGAZ , ZPIZAZ, ZTAUAZ, ZUD  , &
   &  ZFDUV, ZCDUV &
   &)


!-- finally multiply by mu_zero*Incident solar radiation (in relevant spectral interval)

  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFD(JL,JPUV,JK)=ZFDUV(JL,JK)*ZFACT(JL)
      ZCD(JL,JPUV,JK)=ZCDUV(JL,JK)*ZFACT(JL)
     ENDDO
  ENDDO


ELSEIF (KMODTS > 0) THEN
!simplified version of RRTM-SW for monochromatic UV computations
!see in reftra the other possible values of kmodts
!If ZCLEAR=0 the cloud overlap is simply treated passing only PCLDF
!If ZCLEAR is computed with one of the methods above, in SPCVRT the total
!cloud cover is taken into account in an approximate manner weighting the 
!final fluxes

  DO JL=KIDIA,KFDIA
   CALL SPCVRT_SW &
   &         (YDERAD, KLEV, 1, 1, 0, 1, &
   &            PAP(JL,:), PT(JL,:), PAPH(JL,:), ZTH(JL,:), PTL(JL), ZALBD(JL), ZALBP(JL), &
   &            PCLDF(JL,:), ZCLEAR(JL), ZTAU(JL,:), ZCG(JL,:), ZOMEGA(JL,:), ZTAU(JL,:), &
   &            ZTRAY(JL,:),ZTAUAK(JL,:), ZCGAK(JL,:), ZPIZAK(JL,:), ZRMU(JL), ZTAUO3(JL,:), ZFACT(JL), &
   &            ZBBFD(JL,:,JPUV), ZBBFU(JL,:,JPUV), ZBBCD(JL,:,JPUV), ZBBCU(JL,:,JPUV), &
   &            ZBBFDDIR(JL,:,JPUV), ZBBCDDIR(JL,:,JPUV))
  ENDDO

ENDIF



!*         4.1    FILL THE DIAGNOSTIC ARRAYS 
!                 --------------------------

  DO JL = KIDIA,KFDIA
    IF (KMODTS == 0) THEN
      PUVT(JL,JPUV) = ZFD(JL,JPUV,KLEV+1)
      PUVC(JL,JPUV) = MAX( ZCD(JL,JPUV,KLEV+1), 1.E-18_JPRB)
!these are not defined for the old processor
      ZUVDDT(JL,JPUV) = 1.E-18_JPRB 
      ZUVDDC(JL,JPUV) = 1.E-18_JPRB
    ELSEIF (KMODTS > 0) THEN
      PUVT(JL,JPUV) = ZBBFD(JL,1,JPUV) !ARRAYS COME TOP-DOWN FROM RRTM-SW!!
      PUVC(JL,JPUV) = ZBBCD(JL,1,JPUV)
      ZUVDDT(JL,JPUV) = ZBBFDDIR(JL,1,JPUV)
      ZUVDDC(JL,JPUV) = ZBBCDDIR(JL,1,JPUV)
    ENDIF 
    ZTCO3(JL,JPUV)= ZUDAC(JL,1)
  ENDDO
ENDDO






DO JL=KIDIA,KFDIA
  PUVCTOT(JL)=0.0_JPRB
  PUVTTOT(JL)=0.0_JPRB
  ZUVDDTTOT(JL)=0.0_JPRB
  ZUVDDCTOT(JL)=0.0_JPRB
  ZUVTSI(JL) =0.0_JPRB
ENDDO
DO JUV=1,KUV
  DO JL=KIDIA,KFDIA
    PUVCTOT(JL)=PUVCTOT(JL)+PUVC(JL,JUV)*RCIEAS(JUV)
    PUVTTOT(JL)=PUVTTOT(JL)+PUVT(JL,JUV)*RCIEAS(JUV)
!direct incoming and surface flux convluted with ery function
!useful to interpolate hourly values to local noon
    ZUVTSI(JL) =ZUVTSI(JL)+RSUVB(JUV)*PMU0(JL)*RCIEAS(JUV)
    ZUVDDTTOT(JL)=ZUVDDTTOT(JL)+ZUVDDT(JL,JUV)*RCIEAS(JUV)
    ZUVDDCTOT(JL)=ZUVDDCTOT(JL)+ZUVDDC(JL,JUV)*RCIEAS(JUV)
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
!solar zenith angle
  ZUVPOUT(JL,1)=PMU0(JL)
  IF (PMU0(JL) >= RMUZUV) THEN
!surface erythemal radiation total sky
    ZUVPOUT(JL,2)=PUVTTOT(JL)
!surface erythemal radiation clear sky
    ZUVPOUT(JL,3)=PUVCTOT(JL)
!surface direct erythemal radiation total sky
    ZUVPOUT(JL,4)=ZUVDDTTOT(JL)
!surface direct erythemal radiation clear sky
    ZUVPOUT(JL,5)=ZUVDDCTOT(JL)
!TOA erythemal radiation
    ZUVPOUT(JL,6)=ZUVTSI(JL)
  ENDIF
ENDDO
IF(KUV == 24) THEN
  DO JUV=1,KUV
    DO JL=KIDIA,KFDIA
!surface spectral UV radiation total sky
      IF (JUV+6 <= KLEV) ZUVPOUT(JL,JUV+6    )= PUVT(JL,JUV)
!surface spectral UV radiation clear sky
      IF (JUV+6+KUV <= KLEV) ZUVPOUT(JL,JUV+6+KUV)= PUVC(JL,JUV)
    ENDDO
  ENDDO
ENDIF

!------------------------------------------------------------------------------
!
!*         6.6     SET SOME FIELDS TO ZERO
!                  -----------------------


IF(LUVPOUT .AND. .NOT.LIFSMIN .AND. .NOT. LIFSTRAJ) THEN
  VUVP1%P(KIDIA:KFDIA,:) = ZUVPOUT(KIDIA:KFDIA,:)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UVRADI',1,ZHOOK_HANDLE)
END SUBROUTINE UVRADI

