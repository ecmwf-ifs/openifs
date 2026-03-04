! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_RRTM &
 & ( YDEAERATM,YDML_PHY_AER,YDML_GCONF,KIDIA   , KFDIA   , KLON    , KLEV    , KACTAERO, KSTART  , KSTEP , KCLIM, &
 &   PAPH    , PAEROK  , PAER, PRHCL   , &
 &   PAERTAUL, PAERASYL, PAEROMGL, PAERTAUS, PAERASYS, PAEROMGS, &
 &   PAERTAULJ,PAERTAUSJ)

!*** * AER_RRTM* - COMPUTATION OF AEROSOL OPTICAL PROPERTIES

!**   INTERFACE.
!     ----------
!          *AER_RRTM* IS CALLED FROM *RADLSWR*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2010-03-22
!        JJMorcrette 20120910  climatological stratospheric aerosols
!        JJMorcrette 20131001  using MACC-derived aerosol climatology
!        ABozzo      201411    extra variables to output aerosol optical
!                              properties (LDIAGFORCING=True)


!     PURPOSE.
!     --------
!     Computes the optical properties for the various spectral intervals
!       of the RRTM longwave and shortwave radiation schemes from the 
!       GEMS prognostic aerosols

!-----------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD    , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_AEROSOL_MOD , ONLY : MODEL_PHYSICS_AEROSOL_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE YOEAEROP , ONLY : &
  &   ALFS_SU, ALFS_OM, ALFS_DD, ALFS_SS, ALFS_BC &
  & , ASYS_SU, ASYS_OM, ASYS_DD, ASYS_SS, ASYS_BC &
  & , OMGS_SU, OMGS_OM, OMGS_DD, OMGS_SS, OMGS_BC &
  & , ALFL_SU, ALFL_OM, ALFL_DD, ALFL_SS, ALFL_BC &
  & , ASYL_SU, ASYL_OM, ASYL_DD, ASYL_SS, ASYL_BC &
  & , OMGL_SU, OMGL_OM, OMGL_DD, OMGL_SS, OMGL_BC &
  & , ALF_OM, ALF_SU
USE YOEAERATM, ONLY : TEAERATM, NMAXTAER

USE YOMLUN   , ONLY : NULOUT
USE YOMCT3   , ONLY : NSTEP

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)                  ,INTENT(INOUT):: YDEAERATM
TYPE(MODEL_GENERAL_CONF_TYPE)   ,INTENT(INOUT):: YDML_GCONF
TYPE(MODEL_PHYSICS_AEROSOL_TYPE),INTENT(INOUT):: YDML_PHY_AER
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)  :: KACTAERO  ! Number of radiation active aerosol species
INTEGER(KIND=JPIM),INTENT(IN)  :: KSTART, KSTEP
INTEGER(KIND=JPIM),INTENT(IN)  :: KCLIM

REAL(KIND=JPRB)   ,INTENT(IN)  :: PAPH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PAEROK(KLON,KLEV,KACTAERO)
!old aerosol AOT here used only to pass the background type (1 and 6)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PAER(KLON,6,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PRHCL(KLON,KLEV)

REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERTAUL(KLON,KLEV,16), PAERASYL(KLON,KLEV,16), &
  &                               PAEROMGL(KLON,KLEV,16)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERTAUS(KLON,KLEV,14), PAERASYS(KLON,KLEV,14), &
  &                               PAEROMGS(KLON,KLEV,14),PAERTAULJ(KLON,KLEV,KACTAERO), &
  &                               PAERTAUSJ(KLON,KLEV,KACTAERO)



!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: IAER, IBIN, IEFRH, IIRH, ICAER, ITAER, ITYP, IMAXTAER
INTEGER(KIND=JPIM) :: JAER, JK, JL, JTAB, JWAVL
INTEGER(KIND=JPIM) :: INAER(KACTAERO), IRH(KLON,KLEV), ITYPAER(10)

REAL(KIND=JPRB) :: ZALF(KACTAERO), ZAERMSS(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB) :: ZASY(KACTAERO), ZOMG(KACTAERO), ZDP(KLON,KLEV)
REAL(KIND=JPRB) :: ZAERTAUL      , ZAERTAUS      , ZEPSAER

REAL(KIND=JPRB) :: ZTAU_TROP_BKGR, ZTAU_STRAT_BKGR
REAL(KIND=JPRB) :: ZOMG_TROP_BKGR, ZOMG_STRAT_BKGR
REAL(KIND=JPRB) :: ZASY_TROP_BKGR, ZASY_STRAT_BKGR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_RRTM',0,ZHOOK_HANDLE)
ASSOCIATE(NAERO=>YDML_GCONF%YGFL%NAERO, &
 & RRHTAB=>YDML_PHY_AER%YREAERSNK%RRHTAB, &
 & NTYPAER=>YDEAERATM%NTYPAER, &
 & NVOLOPTP=>YDML_PHY_AER%YREAERVOL%NVOLOPTP, &
 & NSTART=>YDML_GCONF%YRRIP%NSTART, LAERCSTR=>YDEAERATM%LAERCSTR)
ZEPSAER=1.E-20_JPRB

!*  --------------------------------------------------------------------
!
!             1.  PREPARATORY
!                 -----------  

!-- units:
!   ------
! ZALF is the extinction coefficient           m2 g-1
! PAEROK is the aerosol mass mixing ratio      kg kg-1 

!-- if MACC prognostic aerosols, i.e. KCLIM=0, copy NMAXTAER and NTYPAER
IF (KCLIM == 0) THEN
  IMAXTAER=NMAXTAER
  ITYPAER(:)=NTYPAER(:)

!-- if aerosols are MACC-derived aerosol climatology KCLIM=1
ELSEIF (KCLIM == 1) THEN
  IMAXTAER=8
  ITYPAER(1)=3
  ITYPAER(2)=3
  ITYPAER(3)=2
  ITYPAER(4)=2
  ITYPAER(5)=1
  ITYPAER(6)=0
  ITYPAER(7)=0
  ITYPAER(8)=0
ENDIF

ICAER=0
DO JAER=1,IMAXTAER
  IF (ITYPAER(JAER) /= 0) THEN
    ITAER=ITYPAER(JAER)
    DO IAER=1,ITAER
      ICAER=ICAER+1
      INAER(ICAER)=JAER*10+IAER
    ENDDO
  ENDIF
ENDDO

!initializations 
IRH(:,:)=1
PAERTAUL(KIDIA:KFDIA,1:KLEV,1:16)=0._JPRB
PAERTAULJ(KIDIA:KFDIA,1:KLEV,1:12)=0._JPRB
PAEROMGL(KIDIA:KFDIA,1:KLEV,1:16)=0._JPRB
PAERASYL(KIDIA:KFDIA,1:KLEV,1:16)=0._JPRB
PAERTAUS(KIDIA:KFDIA,1:KLEV,1:14)=0._JPRB
PAERTAUSJ(KIDIA:KFDIA,1:KLEV,1:12)=0._JPRB
PAEROMGS(KIDIA:KFDIA,1:KLEV,1:14)=0._JPRB
PAERASYS(KIDIA:KFDIA,1:KLEV,1:14)=0._JPRB
ZALF(:)=0._JPRB
ZOMG(:)=0._JPRB
ZASY(:)=0._JPRB
ZTAU_TROP_BKGR=0._JPRB
ZTAU_STRAT_BKGR=0._JPRB
ZOMG_TROP_BKGR=1._JPRB
ZOMG_STRAT_BKGR=1._JPRB
ZASY_TROP_BKGR=0._JPRB
ZASY_STRAT_BKGR=0._JPRB

!-- inputs PAEROK is aerosol mass mixing ratio in kg kg-1 

!-- define RH index from "clear-sky" (not yet!) relative humidity
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDP(JL,JK) =  PAPH(JL,JK)-PAPH(JL,JK-1)
    DO JTAB=1,12
      IF (PRHCL(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO
!-- the effective relative humidity is the low value (20%) assumed for hydrophobic component of OM
IEFRH=3


!-- mass of aerosols in each layer

DO JAER=1,KACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAERMSS(JL,JK,JAER) = MAX( ZEPSAER, PAEROK(JL,JK,JAER)*ZDP(JL,JK)/RG )
    ENDDO
  ENDDO
ENDDO

! !SELECTIVELY PUT AEROSOL MASS TO ZERO TO TEST 
! !INDIVIDUAL CONTRIBUTIONS
! IF (TEST_AER == 1) THEN
! !ONLY ORGANIC MATTER AND SULPHATE
! ZAERMSS(:,:,1:6) = ZAERMSS(:,:,1:6)*ZEPSAER
! ZAERMSS(:,:,9:10) = ZAERMSS(:,:,9:10)*ZEPSAER
! ELSE IF (TEST_AER ==2) THEN
! !ONLY ORGANIC MATTER, SULPHATE AND SEA SALT
! ZAERMSS(:,:,4:6) = ZAERMSS(:,:,4:6)*ZEPSAER
! ZAERMSS(:,:,9:10) = ZAERMSS(:,:,9:10)*ZEPSAER
! ELSE IF (TEST_AER ==3) THEN
!ONLY ORGANIC MATTER, SULPHATE, SEA SALT AND DUST
!ZAERMSS(:,:,9:10) = ZAERMSS(:,:,9:10)*0.25_JPRB !ZEPSAER
! ELSE IF (TEST_AER ==5) THEN
! !ONLY ORGANIC MATTER and SULPHATE coming from Tegen
! ZAERMSS(:,:,7:8) = ZAERMSS(:,:,7:8)*ZEPSAER
! ZAERMSS(:,:,11) = ZAERMSS(:,:,11)*ZEPSAER
! END IF



      
!*  --------------------------------------------------------------------
!
!             2.  LONGWAVE OPTICAL PROPERTIES
!                 ---------------------------  

DO JWAVL=1,16
  PAERTAUL(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
  PAEROMGL(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
  PAERASYL(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
  DO JAER=1,KACTAERO
    ITYP=INAER(JAER)/10
    IBIN=INAER(JAER)-ITYP*10

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:VS,   8=VS,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:cFA,  2:vSU,  1=VS,
!   N.B.: extinction coefficients are in m2 g-1

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        IIRH=IRH(JL,JK)
!        IF (JAER == 8) IIRH=1        ! the hydrophobic OM gets optical properties for RH > 40%
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALFL_SS(IIRH,JWAVL,IBIN)
          ZOMG(JAER)=OMGL_SS(IIRH,JWAVL,IBIN)
          ZASY(JAER)=ASYL_SS(IIRH,JWAVL,IBIN)
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALFL_DD(IBIN,JWAVL)
          ZOMG(JAER)=OMGL_DD(IBIN,JWAVL)
          ZASY(JAER)=ASYL_DD(IBIN,JWAVL)
        ELSEIF (ITYP == 3) THEN
!-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
        IF (IBIN == 2) IIRH=IEFRH
          ZALF(JAER)=ALFL_OM(IIRH,JWAVL)
          ZOMG(JAER)=OMGL_OM(IIRH,JWAVL)
          ZASY(JAER)=ASYL_OM(IIRH,JWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALFL_BC(JWAVL)
          ZOMG(JAER)=OMGL_BC(JWAVL)
          ZASY(JAER)=ASYL_BC(JWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 9) THEN
          ZALF(JAER)=ALFL_SU(IIRH,JWAVL)
          ZOMG(JAER)=OMGL_SU(IIRH,JWAVL)
          ZASY(JAER)=ASYL_SU(IIRH,JWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZASY(JAER)=0._JPRB
          ENDIF
!-- possibility of choosing what optical properties are given to the the volcanic ash
        ELSEIF (ITYP == 8) THEN
          IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
            IIRH=IEFRH
            ZALF(JAER)=ALFL_SU(IIRH,JWAVL)
            ZOMG(JAER)=OMGL_SU(IIRH,JWAVL)
            ZASY(JAER)=ASYL_SU(IIRH,JWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
            ZALF(JAER)=ALFL_BC(JWAVL)
            ZOMG(JAER)=OMGL_BC(JWAVL)
            ZASY(JAER)=ASYL_BC(JWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
            ZALF(JAER)=ALFL_DD(3,JWAVL)    
            ZOMG(JAER)=OMGL_DD(3,JWAVL)
            ZASY(JAER)=ASYL_DD(3,JWAVL)
          ENDIF
        ELSEIF (ITYP == 8) THEN
          ZALF(JAER)=0._JPRB
          ZOMG(JAER)=0._JPRB
          ZASY(JAER)=0._JPRB
        ENDIF
 
        IF (NSTEP == NSTART .AND. (ZALF(JAER) < 0._JPRB&
  &         .OR. ZAERMSS(JL,JK,JAER) < 0._JPRB) ) THEN

          WRITE(NULOUT,FMT='(" AER_RRTM LW JAER ZALF ",2I4,4E12.5)')&
  &             JAER,IIRH,ZALF(JAER),ZOMG(JAER),ZASY(JAER),ZAERMSS(JL,JK,JAER)
        ENDIF

!- ZAERTAUL (ND = g m-2 * m2 g-1)
        ZAERTAUL              = ZAERMSS(JL,JK,JAER) * 1.E+03_JPRB * ZALF(JAER)
IF (JWAVL==6) THEN
! FOR DIAGNOSTICS
        PAERTAULJ(JL,JK,JAER) = ZAERTAUL
ENDIF
        PAERTAUL(JL,JK,JWAVL) = PAERTAUL(JL,JK,JWAVL)+ ZAERTAUL 
        PAEROMGL(JL,JK,JWAVL) = PAEROMGL(JL,JK,JWAVL)+ ZAERTAUL * ZOMG(JAER)
        PAERASYL(JL,JK,JWAVL) = PAERASYL(JL,JK,JWAVL)+ ZAERTAUL * ZOMG(JAER)*ZASY(JAER)
      ENDDO
    ENDDO

  ENDDO

!add tropospheric and stratospheric background species. 
!Based on organic ar 20% UR and sulphate at 20% UR
!the factor 1e3 cancel during the division/multiplication
!sulphates type and no tropospheric bkg is considered
!conversion AOT->mass is done at 550 nm

  IF (.NOT.LAERCSTR) THEN
  DO JK=1,KLEV
     DO JL=KIDIA,KFDIA
        ZTAU_TROP_BKGR=PAER(JL,1,JK)/ALF_OM(3,9)* ALFL_OM(3,JWAVL)
        ZTAU_STRAT_BKGR=PAER(JL,6,JK)/ALF_SU(3,9)* ALFL_SU(3,JWAVL)
        ZOMG_TROP_BKGR=ZTAU_TROP_BKGR*OMGL_OM(3,JWAVL)
        ZOMG_STRAT_BKGR=ZTAU_STRAT_BKGR*OMGL_SU(3,JWAVL)
        ZASY_TROP_BKGR=ZTAU_TROP_BKGR*OMGL_OM(3,JWAVL)*ASYL_OM(3,JWAVL)
        ZASY_STRAT_BKGR=ZTAU_STRAT_BKGR*OMGL_SU(3,JWAVL)*ASYL_SU(3,JWAVL)

        PAERTAUL(JL,JK,JWAVL) = PAERTAUL(JL,JK,JWAVL)&
        & + ZTAU_TROP_BKGR + ZTAU_STRAT_BKGR
        PAEROMGL(JL,JK,JWAVL) = PAEROMGL(JL,JK,JWAVL)&
        & + ZOMG_TROP_BKGR + ZOMG_STRAT_BKGR
        PAERASYL(JL,JK,JWAVL) = PAERASYL(JL,JK,JWAVL)&
        & + ZASY_TROP_BKGR + ZASY_STRAT_BKGR
     ENDDO
  ENDDO
  ENDIF
!-- total optical thickness in each longwave spectral interval: as sum over all 
!     individual aerosol optical thicknesses
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF( PAERTAUL(JL,JK,JWAVL) > 0._JPRB .AND. PAEROMGL(JL,JK,JWAVL) > 0._JPRB) THEN
        PAERASYL(JL,JK,JWAVL)=PAERASYL(JL,JK,JWAVL)/PAEROMGL(JL,JK,JWAVL)
        PAEROMGL(JL,JK,JWAVL)=PAEROMGL(JL,JK,JWAVL)/PAERTAUL(JL,JK,JWAVL)
      ELSE
        PAERASYL(JL,JK,JWAVL)=0._JPRB
        PAEROMGL(JL,JK,JWAVL)=0._JPRB
      ENDIF
    ENDDO
  ENDDO


ENDDO

!*  --------------------------------------------------------------------
!
!             3.  SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------  

DO JWAVL=1,14
  PAERTAUS(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
  PAEROMGS(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
  PAERASYS(KIDIA:KFDIA,1:KLEV,JWAVL)=0._JPRB
  DO JAER=1,KACTAERO
    ITYP=INAER(JAER)/10
    IBIN=INAER(JAER)-ITYP*10

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:VS,   8=xx,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:vFA,  2:vSU,  1=xx,
!   N.B.: extinction coefficients are in m2 g-1

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        IIRH=IRH(JL,JK)
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALFS_SS(IIRH,JWAVL,IBIN)
          ZOMG(JAER)=OMGS_SS(IIRH,JWAVL,IBIN)
          ZASY(JAER)=ASYS_SS(IIRH,JWAVL,IBIN)
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALFS_DD(IBIN,JWAVL)
          ZOMG(JAER)=OMGS_DD(IBIN,JWAVL)
          ZASY(JAER)=ASYS_DD(IBIN,JWAVL)
        ELSEIF (ITYP == 3) THEN
!-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
          IF (IBIN == 2) IIRH=IEFRH
          ZALF(JAER)=ALFS_OM(IIRH,JWAVL)
          ZOMG(JAER)=OMGS_OM(IIRH,JWAVL)
          ZASY(JAER)=ASYS_OM(IIRH,JWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALFS_BC(JWAVL)
          ZOMG(JAER)=OMGS_BC(JWAVL)
          ZASY(JAER)=ASYS_BC(JWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 9) THEN
          ZALF(JAER)=ALFS_SU(IIRH,JWAVL)
          ZOMG(JAER)=OMGS_SU(IIRH,JWAVL)
          ZASY(JAER)=ASYS_SU(IIRH,JWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZASY(JAER)=0._JPRB
          ENDIF
!-- possibility of choosing what optical properties are given to the the volcanic ash
        ELSEIF (ITYP == 8) THEN
          IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
            IIRH=IEFRH
            ZALF(JAER)=ALFS_SU(IIRH,JWAVL)
            ZOMG(JAER)=OMGS_SU(IIRH,JWAVL)
            ZASY(JAER)=ASYS_SU(IIRH,JWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
            ZALF(JAER)=ALFS_BC(JWAVL)
            ZOMG(JAER)=OMGS_BC(JWAVL)
            ZASY(JAER)=ASYS_BC(JWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
            ZALF(JAER)=ALFS_DD(3,JWAVL)    
            ZOMG(JAER)=OMGS_DD(3,JWAVL)
            ZASY(JAER)=ASYS_DD(3,JWAVL)
          ENDIF
        ELSEIF (ITYP == 8) THEN
          ZALF(JAER)=0._JPRB
          ZOMG(JAER)=0._JPRB
          ZASY(JAER)=0._JPRB
        ENDIF

        IF (NSTEP == NSTART .AND. (ZALF(JAER) < 0._JPRB&
  &         .OR. ZAERMSS(JL,JK,JAER) < 0._JPRB) ) THEN 
          WRITE(NULOUT,FMT='(" AER_RRTM SW JAER ZALF ",2I4,4E12.5)')&
  &             JAER,IIRH,ZALF(JAER),ZOMG(JAER),ZASY(JAER),ZAERMSS(JL,JK,JAER) 
        ENDIF

!- ZAERTAUS (ND = g m-2 * m2 g-1)
        ZAERTAUS              = ZAERMSS(JL,JK,JAER) * 1.E+03_JPRB * ZALF(JAER)
IF (JWAVL==10) THEN
! FOR DIAGNOSTICS
        PAERTAUSJ(JL,JK,JAER) = ZAERTAUS
ENDIF

        PAERTAUS(JL,JK,JWAVL) = PAERTAUS(JL,JK,JWAVL)+ ZAERTAUS 
        PAEROMGS(JL,JK,JWAVL) = PAEROMGS(JL,JK,JWAVL)+ ZAERTAUS * ZOMG(JAER)
        PAERASYS(JL,JK,JWAVL) = PAERASYS(JL,JK,JWAVL)+ ZAERTAUS * ZOMG(JAER)*ZASY(JAER)
      ENDDO
    ENDDO

  ENDDO


!add tropospheric and stratospheric background species. 
!Based on organic ar 20% UR and sulphate at 20% UR
!the factor 1e3 cancel during the division/multiplication
!if LAERCSTR is T the strat aer is added in radlswr to the
!sulphates type and no tropospheric bkg is considered
!conversion AOT->mass is done at 550 nm

  IF (.NOT.LAERCSTR) THEN
  DO JK=1,KLEV
     DO JL=KIDIA,KFDIA
        ZTAU_TROP_BKGR=PAER(JL,1,JK)/ALF_OM(3,9)* ALFS_OM(3,JWAVL)
        ZTAU_STRAT_BKGR=PAER(JL,6,JK)/ALF_SU(3,9)* ALFS_SU(3,JWAVL)
        ZOMG_TROP_BKGR=ZTAU_TROP_BKGR*OMGS_OM(3,JWAVL)
        ZOMG_STRAT_BKGR=ZTAU_STRAT_BKGR*OMGS_SU(3,JWAVL)
        ZASY_TROP_BKGR=ZTAU_TROP_BKGR*OMGS_OM(3,JWAVL)*ASYS_OM(3,JWAVL)
        ZASY_STRAT_BKGR=ZTAU_STRAT_BKGR*OMGS_SU(3,JWAVL)*ASYS_SU(3,JWAVL)

        PAERTAUS(JL,JK,JWAVL) = PAERTAUS(JL,JK,JWAVL)&
        & + ZTAU_TROP_BKGR + ZTAU_STRAT_BKGR
        PAEROMGS(JL,JK,JWAVL) = PAEROMGS(JL,JK,JWAVL)&
        & + ZOMG_TROP_BKGR + ZOMG_STRAT_BKGR
        PAERASYS(JL,JK,JWAVL) = PAERASYS(JL,JK,JWAVL)&
        & + ZASY_TROP_BKGR + ZASY_STRAT_BKGR
     ENDDO
  ENDDO
  ENDIF


!-- total optical thickness in each shortwave spectral interval: as sum over all 
!     individual aerosol optical thicknesses
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF( PAERTAUS(JL,JK,JWAVL) > 0._JPRB .AND. PAEROMGS(JL,JK,JWAVL) > 0._JPRB) THEN
        PAERASYS(JL,JK,JWAVL)=PAERASYS(JL,JK,JWAVL)/PAEROMGS(JL,JK,JWAVL)
        PAEROMGS(JL,JK,JWAVL)=PAEROMGS(JL,JK,JWAVL)/PAERTAUS(JL,JK,JWAVL)
      ELSE
        PAERASYS(JL,JK,JWAVL)=0._JPRB
        PAEROMGS(JL,JK,JWAVL)=0._JPRB
      ENDIF
    ENDDO
  ENDDO

ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_RRTM',1,ZHOOK_HANDLE)
END SUBROUTINE AER_RRTM
