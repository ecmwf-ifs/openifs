! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_BDGTMSS &
 & ( YDEAERATM,YDML_PHY_AER,YDCOMPO,YGFL,KIDIA  , KFDIA   , KLON   , KLEV, KACTAERO, KWHAT, KNWAVL, KTWAVL, KSTART, KSTEP, &
 &   PAPH   , PAEROK  , PCAERO , PDZ     , PO3     , PQ      , PRHCL   , PRHO    , PT , PTH, PAPHIF, &
 &   PODSS   , PODDU  , PODOM   , PODBC   , PODSU   , PODNI, PODAM, PODSOA, PODVFA  , PODVSU  , &
 &   PAEPM1 , PAEPM25 , PAEPM10, &
 &   PAERMSS, PAERMSST, PAERTAU, PAERTAUT, PAERTAUI, PAERTAUB, PAERTAUA, PAERTAUF, PAERAOT, &
 &   PAEROMGI,PAERASYI, PAEREXT, PAEREXTD, PAERABST, PEXTRH, PABSRH, PSIGAER, &
 &   PLISIM,  PLISIA,   PLISIS , PLISIT )

!*** * AER_BDGTMSS* - COMPUTATION OF AEROSOL MASS(ES) AND OPTICAL DEPTH(S)

!**   INTERFACE.
!     ----------
!          *AER_BDGTMSS* IS CALLED FROM *CALLPAR*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-03-30
!     A. Benedetti: 2005-04-27 to compute aerosol optical mass from 
!                              optical depth
!     JJMorcrette : 20060721  16 wavelengths
!     JJMorcrette : 20070124  multi-wavelengths implemented
!     JJMorcrette : 20071101  17 wavelengths + absorption AOD
!     JJMorcrette : 20080407  total omega and g
!     JJMorcrette : 20080508  so2/so4
!     JJMorcrette : 20090311  20 wavelengths (+2 for lidar comparisons)
!     JJMorcrette : 20090706  Lidar simulator
!     JJMorcrette : 20101110  10-m Visibility
!     JJMorcrette : 20120501  extinction profiles
!     JJMorcrette : 20130709  absorption optical thickness and height profiles
!     ABozzo      : 20160421  20 wavelengths (+1 for IASI comparisons at 10 microns)
!     SRémy       : 20160530  review PM10/PM2.5/PM1 formulae
!     SRémy       : 20160916  Add nitrates
!     SRémy       : 20171113  Add distinct SOA

!     PURPOSE.
!     --------
!     Depending on the value of KWHAT, either process:
!     - optical thickness from mass mixing ratio (KWHAT =0)
!     - mass mixing ratio from optical thickness (KWHAT/=0)

!     End products are:
!     - a total optical depth for each aerosol type separately        (PAERTAUT) 
!       and a total aerosol optical depth for all aerosols            (PAERTAUI)
!     - an absorption optical depth for each aerosol type separately  (PAERTAUB) 
!       and an absorption aerosol optical depth for all aerosols      (PAERTAUA)
!       and a total aerosol optical depth for all aerosols < 0.5 um   (PAERTAUF)
!     - a vertically integrated mass of aerosol for each aerosol type (PAERMSST)

!     Depending on the value of KWAVL, the optical depths are produced
!     for different wavelengths (in nm)
!         KWVAL=  340,  355,  380,  400,  440,  469,  500,  532,  555, 
!                 645,  670,  800,  858,  865, 1020, 1064, 1240, 1640, 2130

!-----------------------------------------------------------------------

USE MODEL_PHYSICS_AEROSOL_MOD , ONLY : MODEL_PHYSICS_AEROSOL_TYPE
USE YOM_YGFL    ,ONLY : TYPE_GFLD
USE PARKIND1    ,ONLY : JPIM     ,JPRB
USE YOMHOOK     ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST      ,ONLY : RG
USE YOEAERATM   ,ONLY : TEAERATM
USE YOMCOMPO    ,ONLY : TCOMPO
USE YOEAEROP    ,ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA,&
                 &  ASY_SU, ASY_OM, ASY_DD, ASY_SS, ASY_BC, ASY_NI, ASY_AM, ASY_SOA,&
                 &  OMG_SU, OMG_OM, OMG_DD, OMG_SS, OMG_BC, OMG_NI, OMG_AM, OMG_SOA
USE YOMLUN      ,ONLY : NULOUT
USE YOE_AERODIAG,ONLY : NPAERAOT, JPAERAOT_TOTAL, JPAERAOT_NATURAL, JPAERAOT_ANTHRO

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)    ,INTENT(INOUT):: YDEAERATM
TYPE(MODEL_PHYSICS_AEROSOL_TYPE),INTENT(INOUT):: YDML_PHY_AER
TYPE(TCOMPO)      ,INTENT(INOUT):: YDCOMPO
TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)  :: KACTAERO  ! Number of active aerosol species
INTEGER(KIND=JPIM),INTENT(IN)  :: KWHAT, KSTART, KSTEP
INTEGER(KIND=JPIM),INTENT(IN)  :: KNWAVL, KTWAVL(20)

REAL(KIND=JPRB)   ,INTENT(IN)  :: PAPH(KLON,0:KLEV)         , PT(KLON,KLEV)               , PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PAEROK(KLON,KLEV,KACTAERO), PCAERO(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDZ(KLON,KLEV)            , PO3(KLON,KLEV)              , PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PRHCL(KLON,KLEV)          , PRHO(KLON,KLEV)             , PTH(KLON,KLEV+1)

REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERMSS(KLON,KLEV,KACTAERO), PAERTAU(KLON,KLEV,KACTAERO), PAERAOT(KLON,KLEV,NPAERAOT)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERMSST(KLON,KACTAERO)    , PAEREXT(KLON,KLEV,20)      , PAEREXTD(KLON,KLEV,20)   
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERTAUI(KLON,20), PAEROMGI(KLON,20), PAERASYI(KLON,20), PAERTAUT(KLON,KACTAERO,20)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERTAUA(KLON,20), PAERTAUF(KLON,20), PAERTAUB(KLON,KACTAERO,20)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSIGAER(KLON,6)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PLISIM(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PLISIA(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PLISIS(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PLISIT(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PODSS(KLON), PODDU(KLON), PODOM(KLON),PODSOA(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PODBC(KLON), PODSU(KLON), PODNI(KLON), PODAM(KLON), PODVFA(KLON),PODVSU(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAEPM1(KLON),PAEPM25(KLON),PAEPM10(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PEXTRH(KLON,3,3), PABSRH(KLON,3,3), PAERABST(KLON,KLEV,20)

INTEGER(KIND=JPIM), PARAMETER :: IDIAGRH(3) = (/ 1, 5, 6 /) ! 0%, 40%, 50%

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: IBIN, IEFRH, IIRH, INK, ITYP, IWAVL
INTEGER(KIND=JPIM) :: JAER, JK, JL, JTAB, JWAVL, JRH
INTEGER(KIND=JPIM) :: IRH(KLON,KLEV)
LOGICAL :: LLPRINT

REAL(KIND=JPRB) :: ZALF(KACTAERO), ZAERMSS, ZASY(KACTAERO), ZOMG(KACTAERO), ZETA(KLON,KLEV)
REAL(KIND=JPRB) :: ZAEROMGT(KLON,KACTAERO,20)  , ZAERASYT(KLON,KACTAERO,20)

REAL(KIND=JPRB) :: ZAP(KLON,KLEV), ZCLD(KLON,KLEV), ZCO2(KLON,KLEV)
REAL(KIND=JPRB) :: ZDP(KLON,KLEV), ZNO2(KLON,KLEV)
REAL(KIND=JPRB) :: ZEPSIL, ZDZCUM(KLON), ZVISAER(KLON) , ZCNTAB, ZDENSVIS
REAL(KIND=JPRB) :: ZFAC


#include "aer_lidsim.intfb.h"

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_BDGTMSS',0,ZHOOK_HANDLE)
ASSOCIATE(LAERLISI=>YGFL%LAERLISI, LAERVOL=>YDEAERATM%LAERVOL, &
 & LAERNITRATE=>YDCOMPO%LAERNITRATE, NVISWL=>YDEAERATM%NVISWL, &
 & RSS_DRY_MASSFAC=>YDEAERATM%RSS_DRY_MASSFAC, &
 & LAERDUST_NEWBIN=>YDEAERATM%LAERDUST_NEWBIN, &
 & RSS_RH80_MASSFAC=>YDEAERATM%RSS_RH80_MASSFAC, &
 & NWLID=>YDML_PHY_AER%YREAERLID%NWLID, &
 & RRHTAB=>YDML_PHY_AER%YREAERSNK%RRHTAB, RSSGROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%RSSGROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YDML_PHY_AER%YREAERSNK%RSSDENS_RHTAB, &
 & LAERSOA=>YDCOMPO%LAERSOA, &
 & NTYPAER=>YDEAERATM%NTYPAER, &
 & YAERO_DESC=>YDEAERATM%YAERO_DESC, &
 & NVOLOPTP=>YDML_PHY_AER%YREAERVOL%NVOLOPTP, RMMD_SS=>YDML_PHY_AER%YREAERSNK%RMMD_SS)

ZEPSIL=1.E-12_JPRB

LLPRINT=.FALSE.
IF (LLPRINT .AND. KSTEP == KSTART) WRITE(NULOUT,FMT='(" NVOLOPTP=",I3)') NVOLOPTP

!-- units:
!   ------
! ZALF is the extinction coefficient           m2 g-1
! PAEROK is the aerosol mass mixing ratio      kg kg-1 
! PCAERO is the aerosol density                kg cm-3 (KWHAT=0)
! PCAERO is the optical thickness              ND     (KWHAT=2)

IRH(:,:)=1
PAERTAUI(KIDIA:KFDIA,:)=0._JPRB
PAEROMGI(KIDIA:KFDIA,:)=0._JPRB
PAERASYI(KIDIA:KFDIA,:)=0._JPRB
PAERTAUA(KIDIA:KFDIA,:)=0._JPRB
PAERTAUF(KIDIA:KFDIA,:)=0._JPRB
PAERMSS(KIDIA:KFDIA,1:KLEV,1:KACTAERO)= 0.0_JPRB
PAERAOT(KIDIA:KFDIA,1:KLEV,1:NPAERAOT)=0._JPRB
PAEPM1(KIDIA:KFDIA) =0._JPRB
PAEPM25(KIDIA:KFDIA)=0._JPRB
PAEPM10(KIDIA:KFDIA)=0._JPRB
PEXTRH(KIDIA:KFDIA,:,:)=0._JPRB
PABSRH(KIDIA:KFDIA,:,:)=0._JPRB

!-- the effective relative humidity is the low value (20%) assumed for hydrophobic component of OM
IEFRH=3

IF (KWHAT == 0) THEN 
!-- inputs PCAERO is aerosol density in kg cm-3 

!-- define RH index from "clear-sky" (not yet!) relative humidity
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZETA(JL,JK)= (PAPH(JL,JK-1)+PAPH(JL,JK))/(2._JPRB*PAPH(JL,KLEV))
      ZAP(JL,JK)=(PAPH(JL,JK)+PAPH(JL,JK-1))*0.5_JPRB
      ZDP(JL,JK)= PAPH(JL,JK)-PAPH(JL,JK-1)
      DO JTAB=1,12
        IF (PRHCL(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
          IRH(JL,JK)=JTAB
        ENDIF
      ENDDO
    ENDDO
  ENDDO

!-- mass of aerosols in each layer

  DO JAER=1,KACTAERO
    PAERMSST(KIDIA:KFDIA,JAER)=0._JPRB
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
!-- PAEROK in kg kg-1 ; ZAERMSS in kg m-2
        ZAERMSS = PAEROK(JL,JK,JAER)*(PAPH(JL,JK)-PAPH(JL,JK-1))/RG 
        PAERMSS(JL,JK,JAER)= ZAERMSS
!-PAERMSST in kg m-2 
        PAERMSST(JL,JAER)  = PAERMSST(JL,JAER)+PAERMSS(JL,JK,JAER)
      ENDDO
    ENDDO
  ENDDO

!-- compute visibility
   IWAVL=NVISWL
!  IF (LLPRINT) WRITE(NULOUT,FMT='(''Spectral interval ='',I3)') IWAVL
  DO JL=KIDIA,KFDIA
    PSIGAER(JL,6)=0._JPRB
    ZDZCUM(JL)=0._JPRB 
  ENDDO
  DO JK=1,5
    INK=KLEV+1-JK
    DO JL=KIDIA,KFDIA
      PSIGAER(JL,JK)=0._JPRB
      IIRH=IRH(JL,INK)
      ZDENSVIS=PRHO(JL,INK)
      DO JAER=1,KACTAERO
        ITYP=YAERO_DESC(JAER)%NTYP
        IBIN=YAERO_DESC(JAER)%NBIN
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS( IIRH ,IWAVL,IBIN)
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALF_DD(IBIN,IWAVL)
        ELSEIF (ITYP == 3) THEN
!-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
          IF (IBIN == 2) IIRH=IEFRH
          ZALF(JAER)=ALF_OM( IIRH ,IWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP == 6) THEN
          ZALF(JAER)=ALF_NI(IIRH,IWAVL, IBIN)
        ELSEIF (ITYP == 7) THEN
          ZALF(JAER)=ALF_AM(IIRH,IWAVL)
        ELSEIF (ITYP == 8) THEN
          IF (IBIN < 3) THEN
            ZALF(JAER)=ALF_SOA(IIRH,IWAVL, IBIN)
          ELSE
            ZALF=0._JPRB
          ENDIF
        ELSEIF (ITYP ==9) THEN
          IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
            IIRH=IEFRH
            ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
            ZALF(JAER)=ALF_BC(IWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
            ZALF(JAER)=ALF_DD(3,IWAVL)
          ENDIF
        ENDIF
!- PSIGAER (ND: PAEROK in kg kg-1 ; ZALF in m2 g-1 ; ZDENSVIS in kg m-3 ; ZCNTAB in km-1)
        ZCNTAB = PAEROK(JL,INK,JAER) * ZALF(JAER) * ZDENSVIS * 1.E+06_JPRB
        IF (ITYP == 1) ZCNTAB = ZCNTAB * RSS_RH80_MASSFAC
        PSIGAER(JL,JK)=PSIGAER(JL,JK) + ZCNTAB
!        IF (LLPRINT .AND. JL == KIDIA) THEN
!          WRITE(NULOUT,FMT='(I3,5E13.5)') JAER,ZALF(JAER),PAEROK(JL,KLEV,JAER),ZDENSVIS,ZCNTAB,PSIGAER(JL,JK)
!        ENDIF
      ENDDO
      PSIGAER(JL,6)=PSIGAER(JL,6)+PSIGAER(JL,JK)*PDZ(JL,INK)
      ZDZCUM(JL)   =ZDZCUM(JL)   +PDZ(JL,INK)
!-- formulation from Hinkley (1976) Vm = (3.91/Sigma) * (0.55/Lambda)^1.3
!   is here used to get visibility at 0.55 um, where Sigma is the extinction coefficient in km-1
!   0.012 km-1 is the approximate contribution from Rayleigh (could be improved)
!   NB: no contribution by clouds, rain or snow
      ZVISAER(JL)= 3.91_JPRB / MAX(ZEPSIL, PSIGAER(JL,6))
      IF (LLPRINT .AND. JL == KIDIA) THEN
        WRITE(NULOUT,FMT='(''Visib: '',5E13.5)') PSIGAER(JL,6),ZVISAER(JL)
      ENDIF
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PSIGAER(JL,6)=PSIGAER(JL,6)/ZDZCUM(JL)
  ENDDO 
!-- end computation visibility



!-- compute surface extinction coefficient at (550, 440, 670)nm for RH=(40,50)%
  DO JWAVL=1,3
    IWAVL=KTWAVL(5*JWAVL-4) ! Desired wavelengths are indices 1, 6, 11 in the form with 550nm first
    DO JRH=1,3
      DO JL=KIDIA,KFDIA
        DO JAER=1,KACTAERO
          IIRH=IDIAGRH(JRH)
          ITYP=YAERO_DESC(JAER)%NTYP
          IBIN=YAERO_DESC(JAER)%NBIN
          IF (ITYP == 1) THEN
            ZALF(JAER)=ALF_SS( IIRH ,IWAVL,IBIN)
            ZOMG(JAER)=OMG_SS( IIRH ,IWAVL,IBIN)
          ELSEIF (ITYP == 2) THEN
            ZALF(JAER)=ALF_DD(IBIN,IWAVL)
            ZOMG(JAER)=OMG_DD(IBIN,IWAVL)
          ELSEIF (ITYP == 3) THEN
!-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
            IF (IBIN == 2) IIRH=IEFRH
            ZALF(JAER)=ALF_OM( IIRH ,IWAVL)
            ZOMG(JAER)=OMG_OM( IIRH ,IWAVL)
          ELSEIF (ITYP == 4) THEN
            ZALF(JAER)=ALF_BC(IWAVL)
            ZOMG(JAER)=OMG_BC(IWAVL)
          ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
            ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
            ZOMG(JAER)=OMG_SU( IIRH ,IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
            IF (IBIN == 2) THEN
              ZALF(JAER)=0._JPRB
              ZOMG(JAER)=0._JPRB
            ENDIF
          ELSEIF (ITYP == 6) THEN
            ZALF(JAER)=ALF_NI(IIRH,IWAVL, IBIN)
            ZOMG(JAER)=OMG_NI(IIRH,IWAVL, IBIN)
          ELSEIF (ITYP == 7) THEN
            ZALF(JAER)=ALF_AM(IIRH,IWAVL)
            ZOMG(JAER)=OMG_AM(IIRH,IWAVL)
          ELSEIF (ITYP == 8) THEN
            IF (IBIN < 3) THEN
              ZALF(JAER)=ALF_SOA(IIRH,IWAVL,IBIN)
              ZOMG(JAER)=OMG_SOA(IIRH,IWAVL,IBIN)
            ELSE
              ZALF(JAER)=0._JPRB
              ZOMG(JAER)=0._JPRB
            ENDIF
          ELSEIF (ITYP ==9) THEN
            IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
              IIRH=IEFRH
              ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
              ZOMG(JAER)=OMG_SU( IIRH ,IWAVL)
            ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
                ZALF(JAER)=ALF_BC(IWAVL)
                ZOMG(JAER)=OMG_BC(IWAVL)
            ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
              ZALF(JAER)=ALF_DD(3,IWAVL)
              ZOMG(JAER)=OMG_DD(3,IWAVL)
            ENDIF
          ENDIF
          ZCNTAB = PAEROK(JL,KLEV,JAER) * 1.E+03 * ZALF(JAER) * PRHO(JL,KLEV)
          IF (ITYP == 1) ZCNTAB = ZCNTAB * RSS_RH80_MASSFAC
          PEXTRH(JL,JRH,JWAVL) = PEXTRH(JL,JRH,JWAVL) + ZCNTAB
          PABSRH(JL,JRH,JWAVL) = PABSRH(JL,JRH,JWAVL) + ZCNTAB*(1._JPRB-ZOMG(JAER))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!-- end computation surface extinction coefficients

      
      
  DO JWAVL=1,KNWAVL
    IWAVL=KTWAVL(JWAVL)
    PAEREXT(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
    PAEREXTD(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
    PAERABST(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB

    DO JAER=1,KACTAERO
      PAERTAUT(KIDIA:KFDIA,JAER,JWAVL)=0._JPRB
      ZAEROMGT(KIDIA:KFDIA,JAER,JWAVL)=0._JPRB
      ZAERASYT(KIDIA:KFDIA,JAER,JWAVL)=0._JPRB
      PAERTAUB(KIDIA:KFDIA,JAER,JWAVL)=0._JPRB
      ITYP=YAERO_DESC(JAER)%NTYP
      IBIN=YAERO_DESC(JAER)%NBIN

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:BS,   8=VS,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:FA,   1:SB    1=VS,
!   N.B.: extinction coefficients are in m2 g-1

      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          IIRH=IRH(JL,JK)
          ZFAC = 1.0_JPRB
          IF (ITYP == 1) THEN
            ZALF(JAER)=ALF_SS( IIRH ,IWAVL,IBIN)
            ZOMG(JAER)=OMG_SS( IIRH ,IWAVL,IBIN)
            ZASY(JAER)=ASY_SS( IIRH ,IWAVL,IBIN)
            ZFAC = RSS_RH80_MASSFAC
          ELSEIF (ITYP == 2) THEN
            ZALF(JAER)=ALF_DD(IBIN,IWAVL)
            ZOMG(JAER)=OMG_DD(IBIN,IWAVL)
            ZASY(JAER)=ASY_DD(IBIN,IWAVL)
          ELSEIF (ITYP == 3) THEN
!-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
            IF (IBIN == 2) IIRH=IEFRH
            ZALF(JAER)=ALF_OM( IIRH ,IWAVL)
            ZOMG(JAER)=OMG_OM( IIRH ,IWAVL)
            ZASY(JAER)=ASY_OM( IIRH ,IWAVL)
          ELSEIF (ITYP == 4) THEN
            ZALF(JAER)=ALF_BC(IWAVL)
            ZOMG(JAER)=OMG_BC(IWAVL)
            ZASY(JAER)=ASY_BC(IWAVL)
          ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
            ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
            ZOMG(JAER)=OMG_SU( IIRH ,IWAVL)
            ZASY(JAER)=ASY_SU( IIRH ,IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
            IF (IBIN == 2) THEN
              ZALF(JAER)=0._JPRB
              ZOMG(JAER)=0._JPRB
              ZASY(JAER)=0._JPRB
            ENDIF
          ELSEIF (ITYP == 6) THEN
            ZALF(JAER)=ALF_NI(IIRH ,IWAVL,IBIN)
            ZOMG(JAER)=OMG_NI(IIRH ,IWAVL,IBIN)
            ZASY(JAER)=ASY_NI(IIRH ,IWAVL,IBIN)
          ELSEIF (ITYP == 7) THEN
            ZALF(JAER)=ALF_AM(IIRH ,IWAVL)
            ZOMG(JAER)=OMG_AM(IIRH ,IWAVL)
            ZASY(JAER)=ASY_AM(IIRH ,IWAVL)
          ELSEIF (ITYP == 8) THEN
            IF (IBIN < 3) THEN
              ZALF(JAER)=ALF_SOA(IIRH ,IWAVL,IBIN)
              ZOMG(JAER)=OMG_SOA(IIRH ,IWAVL,IBIN)
              ZASY(JAER)=ASY_SOA(IIRH ,IWAVL,IBIN)
            ELSE
              ZALF(JAER)=0._JPRB
              ZOMG(JAER)=0._JPRB
              ZASY(JAER)=0._JPRB
            ENDIF
          ELSEIF (ITYP ==9) THEN
            IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
              IIRH=IEFRH
              ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
              ZOMG(JAER)=OMG_SU( IIRH ,IWAVL)
              ZASY(JAER)=ASY_SU( IIRH ,IWAVL)
            ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
              ZALF(JAER)=ALF_BC(IWAVL)
              ZOMG(JAER)=OMG_BC(IWAVL)
              ZASY(JAER)=ASY_BC(IWAVL)
            ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
              ZALF(JAER)=ALF_DD(3,IWAVL)
              ZOMG(JAER)=OMG_DD(3,IWAVL)
              ZASY(JAER)=ASY_DD(3,IWAVL)
            ENDIF
          ENDIF


!- PAERTAU (ND = g m-2 * m2 g-1), PAERTAUT (ND)
!- For SS, vary depending on RH
          PAERTAU(JL,JK,JAER)  = PAERMSS(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER)
          PAEREXT(JL,JK,JWAVL) = PAEREXT(JL,JK,JWAVL) &
          &                      + PAEROK(JL,JK,JAER) * ZFAC * 1.E+03 * ZALF(JAER) * PRHO(JL,JK)
          PAERABST(JL,JK,JWAVL) = PAERABST(JL,JK,JWAVL) &
          &                      + PAEROK(JL,JK,JAER) * ZFAC * 1.E+03 * ZALF(JAER)*(1._JPRB-ZOMG(JAER)) * PRHO(JL,JK)
          IF (ITYP == 2) THEN
            PAEREXTD(JL,JK,JWAVL) = PAEREXTD(JL,JK,JWAVL) + PAEROK(JL,JK,JAER) * 1.E+03 * ZALF(JAER) * PRHO(JL,JK)
          ENDIF
          PAERTAUT(JL,JAER,JWAVL) = PAERTAUT(JL,JAER,JWAVL)+PAERTAU(JL,JK,JAER)
          ZAEROMGT(JL,JAER,JWAVL) = ZAEROMGT(JL,JAER,JWAVL)+PAERTAU(JL,JK,JAER)*ZOMG(JAER)
          ZAERASYT(JL,JAER,JWAVL) = ZAERASYT(JL,JAER,JWAVL)+PAERTAU(JL,JK,JAER)*ZOMG(JAER)*ZASY(JAER)
!- PAERTAUB is the absorption optical thickness 
          PAERTAUB(JL,JAER,JWAVL) = PAERTAUB(JL,JAER,JWAVL)+PAERTAU(JL,JK,JAER)*(1._JPRB-ZOMG(JAER))
        ENDDO
      ENDDO

!-- total optical thickness: as sum over all individual aerosol optical thicknesses

      DO JL=KIDIA,KFDIA
        PAERTAUI(JL,JWAVL)=PAERTAUI(JL,JWAVL)+PAERTAUT(JL,JAER,JWAVL)
        PAERTAUA(JL,JWAVL)=PAERTAUA(JL,JWAVL)+PAERTAUB(JL,JAER,JWAVL)
        PAEROMGI(JL,JWAVL)=PAEROMGI(JL,JWAVL)+ZAEROMGT(JL,JAER,JWAVL)
        PAERASYI(JL,JWAVL)=PAERASYI(JL,JWAVL)+ZAERASYT(JL,JAER,JWAVL)
      ENDDO
      IF ( (ITYP == 1 .AND. IBIN == 1) .OR. (ITYP == 2 .AND. IBIN == 1) .OR. ITYP >= 3 ) THEN
        DO JL=KIDIA,KFDIA
          PAERTAUF(JL,JWAVL)=PAERTAUF(JL,JWAVL)+PAERTAUT(JL,JAER,JWAVL)
        ENDDO
      ENDIF

!-- for 532 nm, store the profile of optical thickness for the sum of aerosols, 
!   and separately for the natural (SS+DU) and anthropogenic (OM+BC+SU) aerosols
      IF (JWAVL == 9) THEN
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            PAERAOT(JL,JK,JPAERAOT_TOTAL)=PAERAOT(JL,JK,JPAERAOT_TOTAL)+PAERTAU(JL,JK,JAER)
          ENDDO
        ENDDO

!-- ITYP = 1: SS   ; = 2: DU
        IF (ITYP < 3) THEN
          DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
              PAERAOT(JL,JK,JPAERAOT_NATURAL)=PAERAOT(JL,JK,JPAERAOT_NATURAL)+PAERTAU(JL,JK,JAER)
            ENDDO
          ENDDO
        ELSE
!-- ITYP = 3: OM   ; = 4: BC   ; = 5: SU
          DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
              PAERAOT(JL,JK,JPAERAOT_ANTHRO)=PAERAOT(JL,JK,JPAERAOT_ANTHRO)+PAERTAU(JL,JK,JAER)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

    ENDDO

    DO JL=KIDIA,KFDIA
      IF( PAERTAUI(JL,JWAVL) > 0._JPRB .AND. PAEROMGI(JL,JWAVL) > 0._JPRB) THEN
        PAERASYI(JL,JWAVL)=PAERASYI(JL,JWAVL)/PAEROMGI(JL,JWAVL)
        PAEROMGI(JL,JWAVL)=PAEROMGI(JL,JWAVL)/PAERTAUI(JL,JWAVL)
      ELSE
        PAERASYI(JL,JWAVL)=0._JPRB
        PAEROMGI(JL,JWAVL)=0._JPRB
      ENDIF
    ENDDO

  ENDDO


ENDIF

DO JL=KIDIA,KFDIA

!-- store the total optical depths NB: this assumes the default GEMS/MACC aerosol configuration
!-- store the particulate matters for radius <= 1 um, 2.5 um and 10 um

  PAEPM1(JL) = 0._JPRB
  PAEPM25(JL) = 0._JPRB
  PAEPM10(JL) = 0._JPRB

  JAER=0

  IF (NTYPAER(1) == 3) THEN ! Sea-salt (3-bin)
    PODSS(JL) = PAERTAUT(JL,JAER+1,1)+PAERTAUT(JL,JAER+2,1)+PAERTAUT(JL,JAER+3,1)
    PAEPM1(JL) = PAEPM1(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1)*RSS_DRY_MASSFAC
    PAEPM25(JL) = PAEPM25(JL) + (1.00_JPRB*PAEROK(JL,KLEV,JAER+1) + 0.50_JPRB*PAEROK(JL,KLEV,JAER+2))*RSS_DRY_MASSFAC
    PAEPM10(JL) = PAEPM10(JL) + (1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2)+0.05_JPRB*PAEROK(JL,KLEV,JAER+3))) &
                              & *RSS_DRY_MASSFAC
  ELSE
    PODSS(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(1)

  IF (NTYPAER(2) == 3) THEN ! Desert dust (3-bin)
    PODDU(JL) = PAERTAUT(JL,JAER+1,1)+PAERTAUT(JL,JAER+2,1)+PAERTAUT(JL,JAER+3,1)
    IF (LAERDUST_NEWBIN) THEN
      PAEPM1(JL) = PAEPM1(JL) + 0.50_JPRB*PAEROK(JL,KLEV,JAER+1)
      PAEPM25(JL) = PAEPM25(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1) + 0.15_JPRB*PAEROK(JL,KLEV,JAER+2)
      PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2)) + 0.40_JPRB*PAEROK(JL,KLEV,JAER+3)
    ELSE
      PAEPM1(JL) = PAEPM1(JL) + 0.97_JPRB*PAEROK(JL,KLEV,JAER+1) ! SHOULD THIS NOW BE 1.0 ??
      PAEPM25(JL) = PAEPM25(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
      PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2)) + 0.40_JPRB*PAEROK(JL,KLEV,JAER+3)
    ENDIF
  ELSE
    PODDU(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(2)

  IF (NTYPAER(3) == 2) THEN ! Organic matter (hydrophilic & hydrophobic)
    PODOM(JL) = PAERTAUT(JL,JAER+1,1)+PAERTAUT(JL,JAER+2,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
  ELSE
    PODOM(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(3)

  IF (NTYPAER(4) == 2) THEN ! Black carbon (hydrophilic & hydrophobic)
    PODBC(JL) = PAERTAUT(JL,JAER+1,1)+PAERTAUT(JL,JAER+2,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
  ELSE
    PODBC(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(4)

  IF (NTYPAER(5) == 1 .OR. NTYPAER(5) == 2) THEN ! Sulphate (SO4 and possibly SO2 precursor, ignored)
    PODSU(JL) = PAERTAUT(JL,JAER+1,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1)
  ELSE
    PODSU(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(5)

  IF (NTYPAER(6) == 2) THEN ! Nitrate (fine and coarse)
    PODNI(JL) = PAERTAUT(JL,JAER+1,1)+PAERTAUT(JL,JAER+2,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*PAEROK(JL,KLEV,JAER+1) + 0.25_JPRB*PAEROK(JL,KLEV,JAER+2)
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
  ELSE
    PODNI(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(6)

  IF (NTYPAER(7) == 1) THEN ! Ammonium (one bin)
    PODAM(JL) = PAERTAUT(JL,JAER+1,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1)
  ELSE
    PODNI(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(7)

  IF (NTYPAER(8) >=2) THEN ! Secondary organics (2 bins)
    PODSOA(JL) = PAERTAUT(JL,JAER+1,1)+PAERTAUT(JL,JAER+2,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*(PAEROK(JL,KLEV,JAER+1)+PAEROK(JL,KLEV,JAER+2))
  ENDIF
  JAER = JAER + NTYPAER(8)

  IF (NTYPAER(9) == 1) THEN ! Volcanic ash (one bin, treated like consistently with optical properties)
    PODVFA(JL) = PAERTAUT(JL,JAER+1,1)
    IF (NVOLOPTP == 1) THEN ! like sulphate
      PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*PAEROK(JL,KLEV,JAER+1)
      PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*PAEROK(JL,KLEV,JAER+1)
      PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1)
    ELSEIF (NVOLOPTP == 2) THEN ! like black carbon
      PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*PAEROK(JL,KLEV,JAER+1)
      PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*PAEROK(JL,KLEV,JAER+1)
      PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1)
    ELSEIF (NVOLOPTP == 3) THEN ! like coarse dust
      IF (LAERDUST_NEWBIN) THEN
        PAEPM10(JL) = PAEPM10(JL) + 0.40_JPRB*PAEROK(JL,KLEV,JAER+1)
      ELSE
        PAEPM10(JL) = PAEPM10(JL) + 0.14_JPRB*PAEROK(JL,KLEV,JAER+1)
      ENDIF
    ENDIF
  ELSE
    PODVFA(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(9)

  IF (NTYPAER(10) == 1 .OR. NTYPAER(10) == 2) THEN ! Volcanic sulphate (SO4 and possibly SO2 precursor, ignored)
    PODVSU(JL) = PAERTAUT(JL,JAER+1,1)
    PAEPM1(JL) = PAEPM1(JL) + 0.60_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM25(JL) = PAEPM25(JL) + 0.70_JPRB*PAEROK(JL,KLEV,JAER+1)
    PAEPM10(JL) = PAEPM10(JL) + 1.00_JPRB*PAEROK(JL,KLEV,JAER+1)
  ELSE
    PODVSU(JL) = 0._JPRB
  ENDIF
  JAER = JAER + NTYPAER(10)

  ! Convert PM mass mixing ratio to mass concentration
  PAEPM1(JL) = PAEPM1(JL)*PRHO(JL,KLEV)
  PAEPM25(JL) = PAEPM25(JL)*PRHO(JL,KLEV)
  PAEPM10(JL) = PAEPM10(JL)*PRHO(JL,KLEV)
ENDDO

! Watch out the optical properties for this inversion below from optical thickness to mmr might not be 
! correct anymore

! BEN
IF (KWHAT == 2) THEN 
!-- inputs PCAERO are optical thicknesses at 0.55 microns. From those
! get aerosol mass in kg/kg

!-- define RH index from "clear-sky" (not yet!) relative humidity
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      DO JTAB=1,12
        IF (PRHCL(JL,JK)*100. > RRHTAB(JTAB)) THEN
          IRH(JL,JK)=JTAB
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  PAERTAUI(KIDIA:KFDIA,:)=0._JPRB
  PAEROMGI(KIDIA:KFDIA,:)=0._JPRB
  PAERASYI(KIDIA:KFDIA,:)=0._JPRB
  PAERTAUA(KIDIA:KFDIA,:)=0._JPRB
  IWAVL=KTWAVL(1)

  DO JAER=1,KACTAERO
    PAERMSST(KIDIA:KFDIA,JAER)=0._JPRB
    PAERTAUT(KIDIA:KFDIA,JAER,:)=0._JPRB
    PAERTAUB(KIDIA:KFDIA,JAER,:)=0._JPRB

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        
        ITYP=YAERO_DESC(JAER)%NTYP
        IBIN=YAERO_DESC(JAER)%NBIN
!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM,   4:BC,   5:SU,   6:VOL     
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   1:OM
        IIRH=IRH(JL,JK)
        ZFAC = 1.0_JPRB
        IF (ITYP == 3 .AND. IBIN == 2) IIRH=IEFRH
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS( IIRH ,IWAVL,IBIN)
          ZFAC = RSS_RH80_MASSFAC
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALF_DD(IBIN,IWAVL)
        ELSEIF (ITYP == 3) THEN
          ZALF(JAER)=ALF_OM( IIRH ,IWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP == 6) THEN
          ZALF(JAER)=ALF_NI(IIRH ,IWAVL,IBIN)
        ELSEIF (ITYP == 7) THEN
          ZALF(JAER)=ALF_AM(IIRH ,IWAVL)
        ELSEIF (ITYP == 8) THEN
          IF (IBIN < 3) THEN
            ZALF(JAER)=ALF_SOA(IIRH ,IWAVL,IBIN)
          ELSE
            ZALF(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP ==9) THEN
          IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
            IIRH=IEFRH
            ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
            ZALF(JAER)=ALF_BC(IWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
            ZALF(JAER)=ALF_DD(3,IWAVL)
          ENDIF
        ELSEIF (ITYP == 8) THEN
          ZALF(JAER)=0._JPRB
        ENDIF

        PAERMSS(JL,JK,JAER)=0._JPRB

        PAERTAU(JL,JK,JAER)=PCAERO(JL,JK,JAER)
        PAERTAUT(JL,JAER,1)=PAERTAUT(JL,JAER,1)+PAERTAU(JL,JK,JAER)
        
        IF (ZALF(JAER) /= 0.0_JPRB .OR. ZDP(JL,JK) /=0.0_JPRB ) THEN 
          PAERMSS(JL,JK,JAER)=PAERTAU(JL,JK,JAER)*RG/(ZDP(JL,JK)*ZFAC*ZALF(JAER)*1000._JPRB)
        ENDIF

        PAERMSST(JL,JAER)  =PAERMSST(JL,JAER)+PAERMSS(JL,JK,JAER)
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      PAERTAUI(JL,1)=PAERTAUI(JL,1)+PAERTAUT(JL,JAER,1)
      PAERTAUA(JL,1)=0._JPRB
    ENDDO
  ENDDO
ENDIF

!-----------------------------------------------------------------------

!*         2.     LIDAR SIGNAL SIMULATIONS AT 355, 532 and 1064 nm
!                 ------------------------------------------------

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    PLISIM(JL,1:NWLID,JK)=0._JPRB
    PLISIA(JL,1:NWLID,JK)=0._JPRB
    PLISIS(JL,1:NWLID,JK)=0._JPRB
    PLISIT(JL,1:NWLID,JK)=0._JPRB
  ENDDO
ENDDO

IF (LAERLISI) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZCO2(JL,JK)=5.014E-06_JPRB
      ZCLD(JL,JK)=0._JPRB
      ZNO2(JL,JK)=1.E-07_JPRB
    ENDDO
  ENDDO

  CALL AER_LIDSIM &
    &( YDEAERATM,YDML_PHY_AER,YDCOMPO,YGFL, &
    &  KIDIA , KFDIA, KLON, KLEV, KACTAERO, &
    &  PAEROK, PAPH , ZAP , ZCO2, ZCLD    , &
    &  ZDP   , PDZ  , ZNO2, PO3 , PRHCL   , PT , &
    &  PLISIM, PLISIA, PLISIS, PLISIT &
    &)     

ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_BDGTMSS',1,ZHOOK_HANDLE)
END SUBROUTINE AER_BDGTMSS

