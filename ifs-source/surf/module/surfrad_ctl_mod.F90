MODULE SURFRAD_CTL_MOD
CONTAINS
SUBROUTINE SURFRAD_CTL(KDDN,KMMN,KMON,KSECO,&
 & KIDIA,KFDIA,KLON,KTILES,KSW,KLW,&
 & LDNH,&
 & PALBF,PALBICEF,PTVH,&
 & PALCOEFF,PCUR,PCVH,&
 & PASN,PMU0,PTS,PWND,&
 & PWS1,KSOTY,PFRTI,PHLICE,PTLICE,&  
 & YDCST,YDLW,YDSW,YDRAD,YDRDI,YDSOIL,YDFLAKE,&
 & YDURB,PALBD,PALBP,PALB,&
 & PSPECTRALEMISS,PEMIT,&
 & PALBTI,PCCNL,PCCNO,&
 & LNEMOLIMALB)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_LW   , ONLY : TLW
USE YOS_SW   , ONLY : TSW
USE YOS_RAD  , ONLY : TRAD
USE YOS_RDI  , ONLY : TRDI
USE YOS_SOIL , ONLY : TSOIL
USE YOS_FLAKE, ONLY : TFLAKE
USE YOS_URB  , ONLY : TURB
USE CANALB_MOD
USE ABORT_SURF_MOD

! (C) Copyright 1991- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURFRAD  - COMPUTES RADIATIVE PROPERTIES OF SURFACE

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SURFRAD* FROM *CALLPAR*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! LDNH   : LOGICAL  : .TRUE. FOR Northern Hemisphere
! PALBF  : REAL     : FIXED BACKGROUND SURFACE SHORTWAVE ALBEDO
! PALBICEF REAL     : FIXED SEA-ICE ALBEDO (FROM COUPLER)
! PTVH   : REAL     : DOMINANT HIGH VEGETATION TYPE
! PCUR   : REAL     : URBAN COVER (CORRECTED)
! PCVH   : REAL     : HIGH VEGETATION
! PALCOEFF : REAL   : MODIS albedo coefficients:
!   For the 4-component scheme, the second dimension indexes
!      1: UV/Vis direct,  2: UV/Vis diffuse,
!      3: Near-IR direct, 4: Near-IR diffuse,
!   For the 6-component scheme, the second dimension indexes
!      1: UV/Vis isotropic,  2: UV/Vis volumetric,  3: UV/Vis geometric
!      4: Near-IR isotropic, 5: Near-IR volumetric, 6: Near-IR geometric
! PASN   : REAL     : ALBEDO OF EXPOSED SNOW (TYPE 5)
! PMU0   : REAL     : COSINE OF SOLAR ZENITH ANGLE
! PTS    : REAL     : SURFACE TEMPERATURE
! PWND   : REAL     : WIND INTENSITY AT LOWEST LEVEL
! PWS1   : REAL     : TOP LAYER SOIL MOISTURE CONTENT
! KSOTY  : INTEGER  : SOIL TYPE                           (1-7)
! PFRTI  : REAL     : TILE FRACTIONS                      (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
! PHLICE : REAL     : LAKE ICE THICKNESS (m) 
! PTLICE : REAL     : LAKE ICE TEMPERATURE (K) 

!     ==== OUTPUTS ===
! PALBD  : REAL     : SURFACE ALBEDO FOR DIFFUSE RADIATION
! PALBP  : REAL     : SURFACE ALBEDO FOR PARALLEL RADIATION
! PALB   : REAL     : AVERAGE SW ALBEDO (DIAGNOSTIC ONLY)
! PSPECTRALEMISS : REAL : SURFACE LONGWAVE SPECTRAL EMISSIVITY
! PEMIT  : REAL     : BROADBAND SURFACE LONGWAVE EMISSIVITY
! PALBTI : REAL     : BROADBAND ALBEDO FOR TILE FRACTIONS
! PCCNL  : REAL     : CCN CONCENTRATION OVER LAND
! PCCNO  : REAL     : CCN CONCENTRATION OVER OCEAN

!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

! WARNING: this routine is used in several configurations.  In the
! full nonlinear model it is configured with YDRAD%NSW=6 shortwave
! spectral intervals (see susrad_mod.F90), and called requesting
! shortwave albedos in KSW=6 intervals.  In the offline surface scheme
! it is configured with 2 shortwave intervals and called with 2.
! However, in data assimilation, it is configured with 6 intervals,
! but after the nonlinear model has run, this routine is called again
! (by callpartl.F90 and callparad.F90) requesting only KSW=2
! intervals. This results in the albedo returned in the near infrared
! being incorrect, as it is taken from interval 2 of the various
! 6-interval albedo parameterizations in this routine.  This doesn't
! actually matter as callpar??.F90 immediately call sualb2si, which
! overwrites them with the correct values. However, the broadband
! albedo of each tile PALBTI (returned by the present subroutine),
! *is* used in by the subsequent TL/AD physics code, and PALBTI is
! strongly underestimated since it is weighted by the first two RSUN
! elements even though RSUN is configured for 6 intervals. This ought
! to be fixed.  Robin Hogan, 8 May 2019.

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCE.
!     ----------
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     J.-J. MORCRETTE  E.C.M.W.F.    91/03/15

!     MODIFICATIONS.
!     --------------
!     J.-J. MORCRETTE  ECMWF 94/11/15  DIRECT/DIFFUSE ALBEDOS
!     J.-J. MORCRETTE 96/06/07  moisture dep. emissiv. / spectral alb.  
!     PJANSSEN/JJMORCRETTE ECMWF     96/11/07  WIND DEPENDENT SEA ALBEDO
!     PViterbo         ECMWF 99/03/03  Albedo for tile fractions
!     J.-J. Morcrette ECMWF 00/10/24 Spectral albedo for all surfaces
!     JJMorcrette     01-10-08  CCNs concentration over ocean
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     JJMorcrette 20060511 MODIS albedo
!     G. Balsamo    ECMWF   08-01-2006 Include Van Genuchten Hydro.
!     G. Balsamo    ECMWF   03-07-2006 Add soil type
!     E. Dutra/G. Balsamo   01-05-2008 Add lake tile
!     E. Dutra              16-11-2009 snow 2099 cleaning / changes albedo of shaded snow 
!     Linus Magnusson       28-09-2010 Sea-ice
!     Robin Hogan   ECMWF   15-01-2019 MODIS albedo 2x3-components 
!     Robin Hogan   ECMWF   26-02-2019 Removed general spectral rescaling (RWEIGHT)
!     Robin Hogan   ECMWF   26-02-2019 Use Moody et al. for snow albedo in 2 spectral bands
!-----------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)  :: KDDN
INTEGER(KIND=JPIM), INTENT(IN)  :: KMMN
INTEGER(KIND=JPIM), INTENT(IN)  :: KMON(:)
INTEGER(KIND=JPIM), INTENT(IN)  :: KSECO
INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)  :: KSW
INTEGER(KIND=JPIM), INTENT(IN)  :: KLW
INTEGER(KIND=JPIM), INTENT(IN)  :: KLON
INTEGER(KIND=JPIM), INTENT(IN)  :: KTILES
INTEGER(KIND=JPIM), INTENT(IN)  :: KSOTY(:)

REAL(KIND=JPRB),    INTENT(IN)  :: PALBF(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PALBICEF(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PTVH(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PCUR(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PCVH(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PALCOEFF(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PASN(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PMU0(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PTS(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PWND(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PWS1(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PHLICE(:) 
REAL(KIND=JPRB),    INTENT(IN)  :: PTLICE(:)
TYPE(TCST),         INTENT(IN)  :: YDCST
TYPE(TLW),          INTENT(IN)  :: YDLW
TYPE(TSW),          INTENT(IN)  :: YDSW
TYPE(TRAD),         INTENT(IN)  :: YDRAD
TYPE(TRDI),         INTENT(IN)  :: YDRDI
TYPE(TSOIL),        INTENT(IN)  :: YDSOIL
TYPE(TFLAKE),       INTENT(IN)  :: YDFLAKE
TYPE(TURB),         INTENT(IN)  :: YDURB
REAL(KIND=JPRB),    INTENT(OUT) :: PALBD(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PALBP(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PALBTI(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PALB(:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSPECTRALEMISS(:,:) ! dimensioned (KLON,KLW)
REAL(KIND=JPRB),    INTENT(OUT) :: PEMIT(:)
REAL(KIND=JPRB),    INTENT(OUT) :: PCCNL(:)
REAL(KIND=JPRB),    INTENT(OUT) :: PCCNO(:)
LOGICAL        ,    INTENT(IN)  :: LNEMOLIMALB

LOGICAL,   INTENT(IN)  :: LDNH(:)

!     -----------------------------------------------------------------
!*       0.2   LOCAL ARRAYS.
!              -------------

! Indices to PALCOEFF for the 4-component albedo scheme
INTEGER(KIND=JPIM), PARAMETER :: IALUVP = 1, IALUVD = 2, IALNIP = 3, IALNID = 4
! Indices to PALCOEFF for the 6-component albedo scheme
INTEGER(KIND=JPIM), PARAMETER :: IALUVI = 1, IALUVV = 2, IALUVG = 3
INTEGER(KIND=JPIM), PARAMETER :: IALNII = 4, IALNIV = 5, IALNIG = 6

REAL(KIND=JPRB) :: ZBSUR(KLON),ZTS(KLON),ZEXPL(KLON),ZEXPO(KLON),ZQ(KLON)

! Tile emissivity values in each spectral interval
REAL(KIND=JPRB) :: ZEMISSTI(YDRAD%NLWEMISS,KTILES)
REAL(KIND=JPRB) :: ZADTI1,ZADTI2,ZADTI3,ZADTI4,ZADTI5,ZADTI6,ZADTI7,ZADTI8
REAL(KIND=JPRB) :: ZAPTI1,ZAPTI2,ZAPTI3,ZAPTI4,ZAPTI5,ZAPTI6,ZAPTI7,ZAPTI8
REAL(KIND=JPRB) :: ZADTI9,ZAPTI9 
REAL(KIND=JPRB) :: ZADTI10,ZAPTI10

INTEGER(KIND=JPIM) :: JL, JS, JNU, JSW, IM1, IM2

REAL(KIND=JPRB) :: ZBSPE, ZPROP, ZALBICE_AR, ZALBICE_AN &
 & , ZALBD, ZALBP, ZLEN, ZW1, ZW2,ZSTAND,ZWCP  

REAL(KIND=JPRB) :: ZLICE(KLON),ZLWAT(KLON)
REAL(KIND=JPRB) :: ZALICE

REAL(KIND=JPRB) :: ZAI_AR_TOT, ZAI_AN_TOT, ZRW_TOT
REAL(KIND=JPRB) :: ZALBSCALE_AN(KLON),ZALBSCALE_AR(KLON)

! Terms for computing the albedo from the 6-component MODIS
! climatology (Schaaf et al., 2002, Remote Sens. Environ., 83,
! 135-148).
! Solar zenith angle (radians), squared, cubed
REAL(KIND=JPRB) :: ZSZA, ZSZA2, ZSZA3
! Solar-zenith-angle-dependent volumetric and geometric albedo
! coefficients
REAL(KIND=JPRB) :: ZALBVOLCOEFF(KLON), ZALBGEOCOEFF(KLON)
! White-sky volumetric and geometric albedo coefficients
REAL(KIND=JPRB), PARAMETER :: ZALBVOLDIFFCOEFF = 0.189184_JPRB
REAL(KIND=JPRB), PARAMETER :: ZALBGEODIFFCOEFF =-1.377622_JPRB

! Fraction of surface solar energy in UV/Vis and Near-IR
REAL(KIND=JPRB) :: ZSUN_UV, ZSUN_NI

! Urban canyon albedo (THIS USES CTESSEL-URBAN AND WILL BE UPGRADED)
REAL(KIND=JPRB) :: PCANALB(KLON)

REAL(KIND=JPRB) :: ZSPECTRALEMISS(YDRAD%NLWEMISS)

! Two values of broadband and spectral snow albedo used to interpolate
! prognostic broadband snow abledo into spectral intervals
REAL(KIND=JPRB) :: ZALB_SNOW_BB(2), ZALB_SNOW_SPEC(2)

INTEGER(KIND=JPIM) :: IHIGH_VEG_TYPE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFRAD_CTL_MOD:SURFRAD_CTL',0,ZHOOK_HANDLE)

ASSOCIATE(RDAY=>YDCST%RDAY, RTT=>YDCST%RTT, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & NSIL=>YDLW%NSIL, TSTAND=>YDLW%TSTAND, XP=>YDLW%XP, &
 & LCCNL=>YDRAD%LCCNL, LCCNO=>YDRAD%LCCNO, NALBEDOSCHEME=>YDRAD%NALBEDOSCHEME, &
 & NEMISSSCHEME=>YDRAD%NEMISSSCHEME, RCCNLND=>YDRAD%RCCNLND ,RCCNSEA=>YDRAD%RCCNSEA, &
 & RALBSEAD=>YDRDI%RALBSEAD, RALB_SNOW_FOREST=>YDRDI%RALB_SNOW_FOREST, REMISS_DESERT=>YDRDI%REMISS_DESERT, &
 & REMISS_LAND=>YDRDI%REMISS_LAND, REMISS_SNOW=>YDRDI%REMISS_SNOW, REMISS_SEA=>YDRDI%REMISS_SEA, &
 & REMISS_WEIGHT=>YDRDI%REMISS_WEIGHT, REMISS_OLD_WEIGHT=>YDRDI%REMISS_OLD_WEIGHT, &
 & REPALB=>YDRDI%REPALB, NLWEMISS=>YDRAD%NLWEMISS, &
 & LESN09=>YDSOIL%LESN09, LEVGEN=>YDSOIL%LEVGEN, RWCAP=>YDSOIL%RWCAP, &
 & RWRR=>YDURB%RWRR, RROOALB=>YDURB%RROOALB,RURBEMIS=>YDURB%RURBEMIS, &
 & RWCAPM=>YDSOIL%RWCAPM, RWPWP=>YDSOIL%RWPWP, RWPWPM=>YDSOIL%RWPWPM, &
 & RALBICE_AN=>YDSW%RALBICE_AN, RALBICE_AR=>YDSW%RALBICE_AR, RSUN=>YDSW%RSUN, NUVVIS=>YDRAD%NUVVIS)

!     ------------------------------------------------------------------
!*         1.     INITIAL CALCULATIONS
!                 --------------------

! See warning above.
!IF (KSW /= YDRAD%NSW) THEN
!  ! SURFRAD was configured in susrad_mod.F90 to have NSW spectral
!  ! intervals, and so we have an inconsistency
!  CALL ABORT_SURF('SURFRAD_CTL: INPUT KSW DOES NOT MATCH STORED NSW')
!ENDIF

! Compute coefficients for 6-component MODIS albedo scheme
!  Schaaf et al., 2002, Remote Sens. Environ., 83, 135-148.
!  Lucht et al., 2000, IEEE Trans.Geosci.Remote Sens., 38, 997-998.
!  Lucht, 1998: J.Geophys.Res., 103, 8763-8778.
!  Lucht, Lewis, 2000, Int.J.Remote Sens., 21, 81-98.
DO JL=KIDIA,KFDIA
  ZSZA = ACOS(PMU0(JL))
  ZSZA2=ZSZA*ZSZA
  ZSZA3=ZSZA*ZSZA2
  ZALBVOLCOEFF(JL)=-0.007574_JPRB -0.070987_JPRB*ZSZA2 +0.307588_JPRB*ZSZA3
  ZALBGEOCOEFF(JL)=-1.284909_JPRB -0.166314_JPRB*ZSZA2 +0.041840_JPRB*ZSZA3
ENDDO

! Solar fraction in two parts of the spectrum
ZSUN_UV = SUM(RSUN(1:NUVVIS))
IF (KSW > NUVVIS) THEN
  ZSUN_NI = SUM(RSUN(NUVVIS+1:KSW))
ELSE
  ZSUN_NI = 0.0_JPRB
ENDIF

ZSTAND=1.0_JPRB/TSTAND

DO JL=KIDIA,KFDIA
  PALB(JL)=0.0_JPRB
  ZBSUR(JL)=0.0_JPRB
  PEMIT(JL)=0.0_JPRB
  ZTS(JL)=(PTS(JL)-TSTAND)*ZSTAND
ENDDO
DO JNU=1,KTILES
  DO JL=KIDIA,KFDIA
    PALBTI(JL,JNU)=0.0_JPRB
  ENDDO
ENDDO

!* FIND LAKE POINTS WITH ICE COVER 
IF (LEFLAKE) THEN
  WHERE (PHLICE(KIDIA:KFDIA) > RH_ICE_MIN_FLK )
    ZLICE(KIDIA:KFDIA)=1._JPRB
    ZLWAT(KIDIA:KFDIA)=0._JPRB
  ELSEWHERE
    ZLICE(KIDIA:KFDIA)=0._JPRB
    ZLWAT(KIDIA:KFDIA)=1._JPRB
  ENDWHERE
ELSE
  ! With lakes off no ice is considered
  ZLICE(KIDIA:KFDIA)=0._JPRB
  ZLWAT(KIDIA:KFDIA)=1._JPRB
ENDIF

!     ------------------------------------------------------------------       
!*         2.     SURFACE LONGWAVE WINDOW EMISSIVITY
!                 ----------------------------------

!*         2.1    DEPENDENCE OF EMISSIVITY ON SOIL MOISTURE
!                 ------------------------------------------------

DO JL=KIDIA,KFDIA
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKES                 10 : URBAN
! WATER
  ZEMISSTI(:,1) = REMISS_SEA(1:NLWEMISS)
! SEA-ICE
  ZEMISSTI(:,2) = REMISS_SNOW(1:NLWEMISS)
! WET SKIN
  ZEMISSTI(:,3) = REMISS_LAND(1:NLWEMISS)
! LAKES  
  ZEMISSTI(:,9) = ZEMISSTI(:,1)*ZLWAT(JL) + ZEMISSTI(:,2)*ZLICE(JL)  

! DRY SNOW-FREE LOW-VEG
  IF(LEVGEN)THEN
    IF (KSOTY(JL)> 0_JPIM) THEN
      JS=KSOTY(JL)
      ZWCP=1.0_JPRB/(RWCAPM(JS)-RWPWPM(JS))
      ZPROP = MAX(0.0_JPRB,MIN(1.0_JPRB,(PWS1(JL)-RWPWPM(JS))*ZWCP ) )
    ELSE
      ZPROP = 1.0_JPRB
    ENDIF
  ELSE
    ZWCP=1.0_JPRB/(RWCAP-RWPWP)
    ZPROP = MAX(0.0_JPRB, MIN(1.0_JPRB, (PWS1(JL)-RWPWP)*ZWCP ) )
  ENDIF
  ZEMISSTI(:,4) = REMISS_LAND  (1:NLWEMISS)*ZPROP & 
       &        + REMISS_DESERT(1:NLWEMISS)*(1.0_JPRB-ZPROP)
! SNOW ON LOW-VEG+BARE-SOIL
  ZEMISSTI(:,5) = REMISS_SNOW(1:NLWEMISS)
! DRY SNOW-FREE HIGH-VEG
  ZEMISSTI(:,6) = ZEMISSTI(:,4)
! SNOW UNDER HIGH-VEG
  ZEMISSTI(:,7) = ZEMISSTI(:,4)
! BARE SOIL
  IF (PALBF(JL) > 0.30_JPRB) THEN
!   desert emissivity
    ZEMISSTI(:,8) = REMISS_DESERT(1:NLWEMISS)
  ELSE
    ZEMISSTI(:,8) = ZEMISSTI(:,4)
  ENDIF

! SUM OVER FIRST 8 TILES
  ZSPECTRALEMISS = SUM(ZEMISSTI(:,1:8) * SPREAD(PFRTI(JL,1:8),1,NLWEMISS), 2)

! URBAN
  IF (KTILES .GT. 9) THEN
    ZEMISSTI(:,10) = RURBEMIS
    ZSPECTRALEMISS = SUM(ZEMISSTI(:,1:8) * SPREAD(PFRTI(JL,1:8),1,NLWEMISS), 2) +&
    & ZEMISSTI(:,10) * SPREAD(PFRTI(JL,10),1,NLWEMISS)
  ENDIF

  IF (LEFLAKE) THEN
    ZSPECTRALEMISS = ZSPECTRALEMISS + PFRTI(JL,9)*ZEMISSTI(:,9)
  ENDIF
  IF (KLW < NLWEMISS) THEN
    ! Emissivity has been requested in two spectral intervals but
    ! computed in more.  This is most likely because we are in the
    ! simplified physics scheme, which uses one spectral interval for
    ! the infrared window and another for everything else. Perform a
    ! weighted average to get the appropriate emissivity in these two
    ! intervals.
    PSPECTRALEMISS(JL,1) = SUM(ZSPECTRALEMISS * REMISS_OLD_WEIGHT(1:NLWEMISS,1))
    PSPECTRALEMISS(JL,2) = SUM(ZSPECTRALEMISS * REMISS_OLD_WEIGHT(1:NLWEMISS,2))
  ELSE
    ! Emissivity has been requested in the same intervals as it was
    ! computed: copy
    PSPECTRALEMISS(JL,1:NLWEMISS) = ZSPECTRALEMISS
  END IF

ENDDO

!*         2.2    AVERAGE EMISSIVITY OVER LONGWAVE SPECTRUM
!                 -----------------------------------------
!
! Note that when SURFRAD is called from CALLPAR, the broadband
! emissivity PEMIT is used to perform approximate longwave updates
! every gridpoint and timestep, so it must be accurate in order that
! the longwave fluxes are reliable
!

IF (KLW == 2) THEN
  DO JNU=1,NSIL
    DO JL=KIDIA,KFDIA
      ZBSPE=XP(1,JNU)+ZTS(JL)*(XP(2,JNU)+ZTS(JL)*(XP(3,JNU)&
           & +ZTS(JL)*(XP(4,JNU)+ZTS(JL)*(XP(5,JNU)+ZTS(JL)*(XP(6,JNU)&
           & )))))  
      ZBSUR(JL)=ZBSUR(JL)+ZBSPE
      IF (JNU == 3.OR. JNU == 4) THEN
        PEMIT(JL)=PEMIT(JL)+ZBSPE*PSPECTRALEMISS(JL,1)
      ELSE
        PEMIT(JL)=PEMIT(JL)+ZBSPE*PSPECTRALEMISS(JL,2)
      ENDIF
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PEMIT(JL)=PEMIT(JL)/ZBSUR(JL)
  ENDDO
ELSE
  ! Broadband longwave emissivity estimated as a simple weighted
  ! average using a Planck function of 15 degC
  DO JL=KIDIA,KFDIA
    PEMIT(JL) = SUM(PSPECTRALEMISS(JL,1:NLWEMISS)*REMISS_WEIGHT(1:NLWEMISS))
  ENDDO
ENDIF

!     ----------------------------------------------------------------          

!*         3.     SURFACE SHORTWAVE ALBEDO
!                 ------------------------

! Time interpolate the monthly values

IF (KDDN >= 15) THEN
  IM1=KMMN
  IM2=1+MOD(KMMN,12)
  ZLEN=KMON(IM1)*RDAY
  ZW1=REAL((KMON(IM1)-KDDN+15)*RDAY-KSECO,JPRB)/ZLEN
  ZW2=1.-ZW1
ELSE
  IM1=1+MOD(KMMN+10,12)
  IM2=KMMN
  ZLEN=KMON(IM1)*RDAY
  ZW1=REAL((15-KDDN)*RDAY-KSECO,JPRB)/ZLEN
  ZW2=1.-ZW1
ENDIF

! Scaling for sea-ice albedo
ZALBSCALE_AR=1._JPRB
ZALBSCALE_AN=1._JPRB
IF (LNEMOLIMALB) THEN
  ZAI_AR_TOT=0._JPRB
  ZAI_AN_TOT=0._JPRB
  ZRW_TOT=0._JPRB
  DO JSW=1,KSW
      ZAI_AR_TOT=ZAI_AR_TOT+ &
        & RSUN(JSW)*(ZW1*RALBICE_AR(IM1,JSW)+ZW2*RALBICE_AR(IM2,JSW))
      ZAI_AN_TOT=ZAI_AN_TOT+ &
        & RSUN(JSW)*(ZW1*RALBICE_AN(IM1,JSW)+ZW2*RALBICE_AN(IM2,JSW))
      ZRW_TOT=ZRW_TOT+RSUN(JSW)
  ENDDO

  DO JL=KIDIA,KFDIA
    IF (PALBICEF(JL)>0.0_JPRB) THEN
      ZALBSCALE_AR(JL)=PALBICEF(JL)*ZRW_TOT/ZAI_AR_TOT
      ZALBSCALE_AN(JL)=PALBICEF(JL)*ZRW_TOT/ZAI_AN_TOT
    ENDIF
  ENDDO
ENDIF

CALL CANALB(KIDIA,KFDIA,YDURB,PMU0,PCANALB)

DO JSW=1,KSW
  DO JL=KIDIA,KFDIA

!-- MODIS ALBEDO, SPECTRAL AND PARALLEL VS. DIFFUSE

    IF (NALBEDOSCHEME == 1) THEN
! 4-component MODIS scheme
      IF (JSW <= NUVVIS) THEN
!- using MODIS UV-Vis albedo
        ZALBP=(MIN(1._JPRB-REPALB, MAX( PALCOEFF(JL,IALUVP), REPALB) ) )
        ZALBD=(MIN(1._JPRB-REPALB, MAX( PALCOEFF(JL,IALUVD), REPALB) ) )
      ELSE
!- using MODIS Near IR albedo
        ZALBP=(MIN(1._JPRB-REPALB, MAX( PALCOEFF(JL,IALNIP), REPALB) ) )
        ZALBD=(MIN(1._JPRB-REPALB, MAX( PALCOEFF(JL,IALNID), REPALB) ) )
      ENDIF

    ELSEIF (NALBEDOSCHEME == 2) THEN
! 6-component MODIS scheme
      IF (JSW <= NUVVIS) THEN
!- using MODIS UV-Vis albedo
        ZALBP  = PALCOEFF(JL,IALUVI) &
             & + PALCOEFF(JL,IALUVV) * ZALBVOLCOEFF(JL) &
             & + PALCOEFF(JL,IALUVG) * ZALBGEOCOEFF(JL)
        ZALBD  = PALCOEFF(JL,IALUVI) &
             & + PALCOEFF(JL,IALUVV) * ZALBVOLDIFFCOEFF &
             & + PALCOEFF(JL,IALUVG) * ZALBGEODIFFCOEFF
      ELSE
!- using MODIS Near IR albedo
        ZALBP  = PALCOEFF(JL,IALNII) &
             & + PALCOEFF(JL,IALNIV) * ZALBVOLCOEFF(JL) &
             & + PALCOEFF(JL,IALNIG) * ZALBGEOCOEFF(JL)
        ZALBD  = PALCOEFF(JL,IALNII) &
             & + PALCOEFF(JL,IALNIV) * ZALBVOLDIFFCOEFF &
             & + PALCOEFF(JL,IALNIG) * ZALBGEODIFFCOEFF
      END IF
      ! Security bounds
      ZALBP=(MIN(1._JPRB-REPALB, MAX( ZALBP, REPALB) ) )
      ZALBD=(MIN(1._JPRB-REPALB, MAX( ZALBD, REPALB) ) )

    ELSEIF (NALBEDOSCHEME == 3) THEN
! 2-component MODIS scheme (diffuse parts of 4-component scheme)
      IF (JSW <= NUVVIS) THEN
!- using MODIS UV-Vis albedo
        ZALBD=(MIN(1._JPRB-REPALB, MAX( PALCOEFF(JL,IALUVD), REPALB) ) )
      ELSE
!- using MODIS Near IR albedo
        ZALBD=(MIN(1._JPRB-REPALB, MAX( PALCOEFF(JL,IALNID), REPALB) ) )
      ENDIF
        ZALBP=ZALBD

    ELSE
!-- ERBE ALBEDO, FLAT, ISOTROPIC
      ZALBD = PALBF(JL)
      ZALBP = PALBF(JL)
    ENDIF

! Copy albedos to different surface types

! WET SKIN
    ZADTI3=ZALBD
    ZAPTI3=ZALBP
! DRY SNOW-FREE LOW-VEG
    ZADTI4=ZALBD
    ZAPTI4=ZALBP
! DRY SNOW-FREE HIGH-VEG
    ZADTI6=ZALBD
    ZAPTI6=ZALBP
! BARE SOIL
    ZADTI8=ZALBD
    ZAPTI8=ZALBP
! URBAN
    IF (KTILES .GT. 9) THEN
      ZADTI10=RROOALB*(1.0_JPRB/(1.0_JPRB+RWRR)) + PCANALB(JL)*(RWRR/(1.0_JPRB+RWRR))
      ZAPTI10=RROOALB*(1.0_JPRB/(1.0_JPRB+RWRR)) + PCANALB(JL)*(RWRR/(1.0_JPRB+RWRR))
      ZADTI3=ZADTI3*(1.0_JPRB-PCUR(JL))+ZADTI10*PCUR(JL)    
      ZAPTI3=ZAPTI3*(1.0_JPRB-PCUR(JL))+ZAPTI10*PCUR(JL)
    ENDIF
! WATER
! Taylor et al. for sea
    ZADTI1=RALBSEAD
!    ZADTI1 = 0.07_JPRB !FOR DEBUGGING PURPOSES ONLY
    ZAPTI1=MAX(0.037_JPRB/(1.1_JPRB*PMU0(JL)**1.4_JPRB+0.15_JPRB),REPALB)
!    ZAPTI1=0.08_JPRB    ! DEBUGGING PURPOSES ONLY
    
! SEA-ICE
    ZALBICE_AR=ZW1*RALBICE_AR(IM1,JSW)+ZW2*RALBICE_AR(IM2,JSW)
    ZALBICE_AN=ZW1*RALBICE_AN(IM1,JSW)+ZW2*RALBICE_AN(IM2,JSW)
    !KW! scale with albedo from coupler
    ZALBICE_AR=MIN(1.0_JPRB,ZALBSCALE_AR(JL)*ZALBICE_AR)
    ZALBICE_AN=MIN(1.0_JPRB,ZALBSCALE_AN(JL)*ZALBICE_AN)
    IF (LDNH(JL)) THEN
      ZADTI2=ZALBICE_AR
      ZAPTI2=ZALBICE_AR
    ELSE
      ZADTI2=ZALBICE_AN
      ZAPTI2=ZALBICE_AN
    ENDIF

! SNOW ON LOW-VEG+BARE-SOIL
    
    ! Moody et al. (RSE 2007) observed the following snow albedos:
    !               0.3-0.7micron  0.7-5.0micron  Broadband
    ! Savanna          0.57           0.39          0.47
    ! Permanent snow   0.89           0.57          0.74
    !
    ! However, the link from spectral to broadband is not exactly
    ! consistent with the "RSUN" spectral split here.  Here use the
    ! prognostic snow albedo, PASN, to interpolate between these
    ! values to obtain the albedos in the UV/Vis, with extrapolation,
    ! but preventing values above 0.98. The Near-IR values are then
    ! calculated to ensure the average is consistent with RSUN.
    ZALB_SNOW_BB   = [0.47_JPRB, 0.74_JPRB] ! Broadband reference values
    ZALB_SNOW_SPEC = [0.57_JPRB, 0.89_JPRB] ! UV/Vis reference values
    ! Linear interpolation to obtain albedo in UV/Vis
    ZADTI5 = (ZALB_SNOW_SPEC(1) * (ZALB_SNOW_BB(2) - PASN(JL)) &
         &   +ZALB_SNOW_SPEC(2) * (PASN(JL) - ZALB_SNOW_BB(1))) &
         &  /(ZALB_SNOW_BB(2) - ZALB_SNOW_BB(1))
    ! Security
    ZADTI5 = MIN(0.98_JPRB, MAX(0.02_JPRB, ZADTI5))
    IF (JSW > NUVVIS) THEN
      ! Select Near-IR value to ensure average value is equal to PASN
      ZADTI5 = (PASN(JL) - ZADTI5*ZSUN_UV) / ZSUN_NI
      ZADTI5 = MIN(0.98_JPRB, MAX(0.02_JPRB, ZADTI5)) ! Security
    ENDIF
    ZAPTI5 = ZADTI5 ! Direct albedo = Diffuse albedo

! SNOW UNDER HIGH-VEG
    IF ( LESN09 ) THEN
      IHIGH_VEG_TYPE = NINT(PTVH(JL))
    ELSE
      ! Note that RALB_SNOW_FOREST(0,:) contains the default
      ! snow/forest albedos
      IHIGH_VEG_TYPE = 0
    ENDIF

    IF (JSW <= NUVVIS) THEN
      ZADTI7=RALB_SNOW_FOREST(IHIGH_VEG_TYPE,1) ! UV/Vis
    ELSE
      ZADTI7=RALB_SNOW_FOREST(IHIGH_VEG_TYPE,2) ! Near-IR
    END IF
    ZAPTI7 = ZADTI7 ! Direct albedo = Diffuse albedo

! CORRECT FOR URBAN - TILE 7 (URBAN SNOW = HIGVEG SNOW)
    IF (KTILES .GT. 9) THEN
     IF (PCUR(JL) .GT. 0.0_JPRB) THEN
      ZAPTI7 = ZAPTI7*(PCVH(JL)*(1.0_JPRB-PCUR(JL))) &
      & +(0.5_JPRB*ZAPTI7+0.5_JPRB*ZAPTI10)*PCUR(JL)
      ZAPTI7 = ZAPTI7/(PCVH(JL)*(1.0_JPRB-PCUR(JL)) + PCUR(JL))
      ZADTI7 = ZADTI7*(PCVH(JL)*(1.0_JPRB-PCUR(JL))) &
      & +(0.5_JPRB*ZADTI7+0.5_JPRB*ZADTI10)*PCUR(JL)
      ZADTI7 = ZADTI7/(PCVH(JL)*(1.0_JPRB-PCUR(JL)) + PCUR(JL))
     ENDIF

    ENDIF

!  LAKES 
    !* use new formulation for lake ice albedo calculation 
    IF (LEFLAKE) THEN
      ZALICE=EXP(-95.6_JPRB*(RTT-PTLICE(JL))/RTT)
    ELSE
      ZALICE=1._JPRB
    ENDIF
    ZALICE=0.7_JPRB*(1.0_JPRB-ZALICE)+0.4_JPRB*ZALICE
    
    !* ZALICE=ZADTI2 !* use montlhy sea-ice albedo 
    
    ZADTI9=ZADTI1*ZLWAT(JL)+ZALICE*ZLICE(JL)
    ZAPTI9=ZAPTI1*ZLWAT(JL)+ZALICE*ZLICE(JL)

    
!* AVERAGE (DIFFUSE+PARALLEL) ALBEDO OVER THE WHOLE SHORTWAVE SPECTRUM, FOR 
! FOR EACH TILE (USED IN SURFACE AND VERT.DIFFUSION SCHEMES).
! See warning above: this returns incorrect results in TL/AD since
! KSW=2 but YRERAD%NSW=6, so the integral over RSUN is incomplete.
    PALBTI(JL,1)=PALBTI(JL,1)+RSUN(JSW)*(ZADTI1+ZAPTI1)*0.5_JPRB
    PALBTI(JL,2)=PALBTI(JL,2)+RSUN(JSW)*(ZADTI2+ZAPTI2)*0.5_JPRB
    PALBTI(JL,3)=PALBTI(JL,3)+RSUN(JSW)*(ZADTI3+ZAPTI3)*0.5_JPRB
    PALBTI(JL,4)=PALBTI(JL,4)+RSUN(JSW)*(ZADTI4+ZAPTI4)*0.5_JPRB
    PALBTI(JL,5)=PALBTI(JL,5)+RSUN(JSW)*(ZADTI5+ZAPTI5)*0.5_JPRB
    PALBTI(JL,6)=PALBTI(JL,6)+RSUN(JSW)*(ZADTI6+ZAPTI6)*0.5_JPRB
    PALBTI(JL,7)=PALBTI(JL,7)+RSUN(JSW)*(ZADTI7+ZAPTI7)*0.5_JPRB
    PALBTI(JL,8)=PALBTI(JL,8)+RSUN(JSW)*(ZADTI8+ZAPTI8)*0.5_JPRB
    IF (KTILES .GT. 9) THEN
     PALBTI(JL,10)=PALBTI(JL,10)+RSUN(JSW)*(ZADTI10+ZAPTI10)*0.5_JPRB
    ENDIF
    IF (LEFLAKE) THEN 
      PALBTI(JL,9)=PALBTI(JL,9)+RSUN(JSW)*(ZADTI9+ZAPTI9)*0.5_JPRB   
    ENDIF
    
!* SUM OVER TILES SPECTRALLY DEFINED FOR THE INDIVIDUAL SPECTRAL 
!  INTERVALS OF THE SHORTWAVE RADIATION SCHEME. WEIGHTING IS BY 
!  FRACTIONAL COVER OF EACH TILE.
!  NB: THE DISTINCTION BETWEEN DIRECT AND DIFFUSE IS ACTUALLY 
!  IRRELEVANT FOR PRE-MODIS ERA ALBEDO   
!
!* DONE FOR DIFFUSE ALBEDO
!
    PALBD(JL,JSW)=PFRTI(JL,1)*ZADTI1 &
     & +PFRTI(JL,2)*ZADTI2 &
     & +PFRTI(JL,3)*ZADTI3 &
     & +PFRTI(JL,4)*ZADTI4 &
     & +PFRTI(JL,5)*ZADTI5 &
     & +PFRTI(JL,6)*ZADTI6 &
     & +PFRTI(JL,7)*ZADTI7 &
     & +PFRTI(JL,8)*ZADTI8 

    IF (KTILES .GT. 9) THEN
     PALBD(JL,JSW)=PALBD(JL,JSW)+PFRTI(JL,10)*ZADTI10
    ENDIF

    IF (LEFLAKE) THEN
      PALBD(JL,JSW)=PALBD(JL,JSW) &
     & +PFRTI(JL,9)*ZADTI9    
    ENDIF

!* DONE FOR DIRECT ALBEDO
!
    PALBP(JL,JSW)=PFRTI(JL,1)*ZAPTI1 &
     & +PFRTI(JL,2)*ZAPTI2 &
     & +PFRTI(JL,3)*ZAPTI3 &
     & +PFRTI(JL,4)*ZAPTI4 &
     & +PFRTI(JL,5)*ZAPTI5 &
     & +PFRTI(JL,6)*ZAPTI6 &
     & +PFRTI(JL,7)*ZAPTI7 &
     & +PFRTI(JL,8)*ZAPTI8

    IF (KTILES .GT. 9) THEN
     PALBP(JL,JSW)=PALBP(JL,JSW)+PFRTI(JL,10)*ZAPTI10
    ENDIF

    IF (LEFLAKE) THEN
      PALBP(JL,JSW)=PALBP(JL,JSW) &
     & +PFRTI(JL,9)*ZAPTI9
    ENDIF
     
!* AVERAGE (DIFFUSE) ALBEDO OVER SHORTWAVE SPECTRUM 
!  DIAGNOSTICS ONLY FOR STAND-ALONE SURFACE SCHEME,  
!  SUPERSEDED IN FULL MODEL BY DIRECT+DIFFUSE ONE IN *CALLPAR*
!
    PALB(JL)=PALB(JL)+RSUN(JSW)*PALBD(JL,JSW)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         4.     CCN EMISSION OVER OCEAN
!                 -----------------------

!* representation of surface concentration from Genthon, 1992, Tellus 4
IF (LCCNO) THEN      
  DO JL=KIDIA,KFDIA
    IF (PWND(JL) > 30._JPRB) THEN
      ZQ(JL)=327._JPRB  
    ELSEIF (PWND(JL) > 15._JPRB) THEN  
      ZQ(JL)=EXP(0.13_JPRB*PWND(JL)+1.89_JPRB)
    ELSE
      ZQ(JL)=EXP(0.16_JPRB*PWND(JL)+1.44_JPRB)
    ENDIF  
    ZEXPO(JL)=1.2_JPRB+0.5_JPRB*LOG10(ZQ(JL))
    PCCNO(JL)=10._JPRB**ZEXPO(JL)
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    PCCNO(JL)=RCCNSEA
  ENDDO
ENDIF    

IF (LCCNL) THEN      
  DO JL=KIDIA,KFDIA
    IF (PWND(JL) <= 15._JPRB) THEN
      ZQ(JL)=EXP(0.16_JPRB*PWND(JL)+1.45_JPRB)
    ELSE
      ZQ(JL)=EXP(0.13_JPRB*PWND(JL)+1.89_JPRB)
    ENDIF
    ZEXPL(JL)=2.21_JPRB+0.3_JPRB*LOG10(ZQ(JL))
    PCCNL(JL)=10._JPRB**ZEXPL(JL)
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    PCCNL(JL)=RCCNLND
  ENDDO
ENDIF  

!DO JL=KIDIA,KFDIA
!  9401 FORMAT(1x,'CCN ',I6,3E13.5)  
!ENDDO
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFRAD_CTL_MOD:SURFRAD_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SURFRAD_CTL
END MODULE SURFRAD_CTL_MOD
