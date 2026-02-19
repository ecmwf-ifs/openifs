! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE SLTEND( YDMODEL,KIDIA, KFDIA, KLON, KLEV, KTEND,&
 & PT, PQ, PA, PCLD, PRSF1, PRS1, &
 & PDYNU, PDYNV, PDYNT, PDYNQ, PDYNO3, PDYNA, PDYNCLD, &
 & PCENU, PCENV, PCENT, PCENQ, PCENO3, PCENA, PCENCLD, &
 & PTENU, PTENV, PTENT, PTENQ, PTENO3, PTENA, PTENCLD, &
 & PPHTENU, PPHTENV, PPHTENT, PPHTENGFL, &
 & PVFU,  PVFV,  PVFT,  PVFQ,  PVFA, &
 & PSATT, PSATQ, PSATA, PSATCLD, &
 & PSAVTEND, PGFLSLP, PEXTRA, KFLDX)  


!**** *SLTEND * - compute sl tendencies for next time step and the 
!                 updates of the current profile. 
!               - remove supersaturation on final profiles and 
!                 save it for the next time step
!     PURPOSE.
!     --------
!               second order physics, computing the total dynamical and physical 
!               tendencies and the contribution of te physical tendencies
!               to be used in the next time step as part of the first guess
!               and as part of the final profile.

!**   Interface.
!     ----------
!        *CALL* *SLTEND*

! INPUT:

! KIDIA      : START OF HORIZONTAL LOOP
! KFDIA      : END   OF HORIZONTAL LOOP
! KLON       : HORIZONTAL DIMENSION
! KLEV       : END OF VERTICAL LOOP AND VERTICAL DIMENSION
! KTEND      : NUMBER OF PARAMETERS IN PSAVTEND
! KTRAC      : NUMBER OF ACTIVE TRACERS IN PVFC: NGHG+NAER+NGFL_EXT

! PT         : TEMPERATURE AT TIME t
! PQ         : HUMIDITY AT TIME t
! PA         : CLOUD FRACTION AT TIME t

! PRSF1      : PROVISIONAL T+DT PRESSURE ON FULL LEVELS
! PRS1       : PROVISIONAL T+DT PRESSURE ON HALF LEVELS

! time level (t+1)-(t)
! PDYNU      : EUL. DYNAMICAL TENDENCY OF U-COMP. OF WIND.
! PDYNV      : EUL. DYNAMICAL TENDENCY OF V-COMP. OF WIND.
! PDYNT      : EUL. DYNAMICAL TENDENCY OF TEMPERATURE.
! PDYNQ      : EUL. DYNAMICAL TENDENCY OF HUMIDITY
! PDYNO3     : EUL. DYNAMICAL TENDENCY OF OZONE MIXING RATIO (ECMWF PROG. OZONE)
! PDYNA      : EUL. DYNAMICAL TENDENCY OF CLOUD FRACTION

! (note: the updated total cloud liquid water and ice are computed in CLOUDSC !)

! time level (t) at departure point
! PPHTENU    : PHYSICAL TENDENCY OF U-COMP. OF WIND
! PPHTENV    : PHYSICAL TENDENCY OF V-COMP. OF WIND
! PPHTENT    : PHYSICAL TENDENCY OF TEMPERATURE
! PPHTENGFL    : PHYSICAL TENDENCY OF GFL FIELDS with LPHY=T

! arriv. point only at time level (t+1)
! PVFU       : ARRIV. POINT TENDENCY FOR VDIF + GWDRAG OF U-COMP. OF WIND
! PVFV       : ARRIV. POINT TENDENCY FOR VDIF + GWDRAG OF V-COMP. OF WIND
! PVFT       : ARRIV. POINT TENDENCY FOR VDIF OF TEMPERATURE
! PVFQ       : ARRIV. POINT TENDENCY FOR VDIF OF HUMIDITY
! PVFA       : ARRIV. POINT TENDENCY FOR VDIF OF CLOUD FRACTION

! arriv. point only at time level (t+1)
! PSATT       : ARRIV. POINT TENDENCY FOR SATURATION ADJUSTMENT FOR TEMPERATURE
! PSATQ       : ARRIV. POINT TENDENCY FOR SATURATION ADJUSTMENT FOR HUMIDITY
! PSATA       : ARRIV. POINT TENDENCY FOR SATURATION ADJUSTMENT FOR CLOUD FRACTION
! PSATCLD     : ARRIV. POINT TENDENCY FOR SATURATION ADJUSTMENT FOR CLOUD VARIABLES

! UPDATED:

! input : arriv. point only at time level (t+1),
! PCENU      : TENDENCY OF U-COMP. OF WIND.
! PCENV      : TENDENCY OF V-COMP. OF WIND.
! PCENT      : TENDENCY OF TEMPERATURE.
! PCENQ      : TENDENCY OF HUMIDITY
! PCENO3     : TENDENCY OF OZONE MIXING RATIO (ECMWF PROG. OZONE)
! PCENA      : TENDENCY OF CLOUD FRACTION

! output: complete arriv. + departure point tendencies including the dynamical tendencies
! PTENU      : TENDENCY OF U-COMP. OF WIND.
! PTENV      : TENDENCY OF V-COMP. OF WIND.
! PTENT      : TENDENCY OF TEMPERATURE.
! PTENQ      : TENDENCY OF HUMIDITY
! PTENO3     : TENDENCY OF OZONE MIXING RATIO (ECMWF PROG. OZONE)
! PTENA      : TENDENCY OF CLOUD FRACTION
! PTENCLV    : TENDENCY OF CLOUD/PRECIP

! PSAVTEND   : ARRAY OF TENDENCIES + AUX. SUPERSATURATION TO BE SAVED FOR NEXT TIME STEP

!-----------------------------------------------------------------------
!     Method. See documentation.
!   
!       This subroutine returns the following difference of slow physics tendencies for a field X:
!              PTENX=0.5*[-(dX^(+)/dt)_slowPhys+(dX^(-)/dt)_slowphys|D] where D denotes interpolation
!                                                                       to the departure point.
!       The above quantity is added in total (cummulative) tendency PCENX when update_state() which is 
!       called immediately after to return the required tendency when LSLPHY=T:
!              PTENX+PCENX=(dX^(+)/dt)_fastPhys + 0.5*[(dX^(+)/dt)_slowPhys + (dX^(-)/dt)_slowphys|D]   
!
!     -------
!     Modifications.
!     --------------
!     ORIGINAL 2001-06-25, Nils Wedi
!     M.Hamrud    : 03-08-01 GFL introduction
!     M.Hamrud    : 01-Oct-2003 CY28 Cleaning
!     A.Untch     : March-2004 Introduced EXTRA GFL fields (GFL_EXT)
!     M.Ko"hler   : 03-12-2004 VDF cloud tendency added for moist
!                              advection-diffusion PBL
!     A.Untch     : 12-03-2005 Aerosols as named GFL fields
!     J.Flemming  : 11-04-2005 Aerosols replaced with reactive gases
!     A.Tompkins  : CY31R2     Changes for supersaturation code
!     P.Bechtold  : 11-12-2005 Reorganize GFL/Tracer part
!     M.Janiskova : 21-12-2005 Modified condition for using cloud tendency
!     D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!     S. Serrar     07-Sep-2006 tracers added for diagnostics (GEMS)
!     P.Bechtold    04-10-2006 Remove GFL tracers as SLTEND does not give positive definit results
!     A.Tompkins/R.Forbes : 15-03-2010  New CLV cloud variables added
!     R.Forbes    : 01-03-2011 Changed supersat liquid temperature threshold
!     R.Forbes    : 01-10-2011 Limited supersaturation to avoid excessive values
!     K. Yessad (July 2014): Move some variables.
!     M.Diamantakis/F. Vana  : 18-10-2013 More freedom to use GFL attributes and LSLPHY independently
!     J.Hague : 21-10-2014 Vector Optimisation for Cray
!     F. Vana  05-Mar-2015  Support for single precision
!     S. Malardel: Nov 2017 Fix pointers for PPHTENGFL to allow LSLPHY if ICI
!     M. Diamantakis (June 2018): Add LPHY control for q
!     R. Forbes 01-Nov-2018  Pass in arrays for cloud budget update
!     R. Forbes 01-Nov-2018  Implement cloud/precip SL-AV supersat consistently  
!     R. Forbes 01-Feb-2020  Rewrite and simplify the code. Consistent options for cloud vars.
!                            Call to cloud_supersatcheck to remove supersaturation here
!     F. Vana   14-Sep-2020 : More flexibility to implicitness & better memory handling
!     R. Forbes 01-Jul-2021  Included subtraction of saturation adjustment tendencies (fast process) 
!-----------------------------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RD
USE YOMCT3    , ONLY : NSTEP
USE YOECLDP   , ONLY : NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLD(KLON,KLEV,NCLV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSF1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRS1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNCLD(KLON,KLEV,NCLV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENCLD(KLON,KLEV,NCLV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENCLD(KLON,KLEV,NCLV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIMSLP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSATT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSATQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSATA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSATCLD(KLON,KLEV,NCLV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSAVTEND(KLON,KLEV,KTEND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLSLP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIMSLP) 
! Extra fields for diagnostics
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX ! Number of extra fields

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL, IU, IV, IT, IS

! Controls T-dependence for liquid/ice production (1=mixedphase, 2=homogfrz)
INTEGER(KIND=JPIM) :: IFTLIQICE 

! required for supersaturation adjustments
REAL(KIND=JPRB) :: ZADJA(KLON,KLEV)
REAL(KIND=JPRB) :: ZADJL(KLON,KLEV)
REAL(KIND=JPRB) :: ZADJI(KLON,KLEV)
REAL(KIND=JPRB) :: ZT_UPD(KLON)     ! Updated temperature for supersat check
REAL(KIND=JPRB) :: ZQ_UPD(KLON)     ! Updated humidity for supersat check
REAL(KIND=JPRB) :: ZA_UPD(KLON)     ! Updated cloud fraction for supersat check
REAL(KIND=JPRB) :: ZT_ADJ(KLON)     ! Supersat check temperature change
REAL(KIND=JPRB) :: ZQ_ADJ(KLON)     ! Supersat check humidity change
REAL(KIND=JPRB) :: ZA_ADJ(KLON)     ! Supersat check cloud fraction change
REAL(KIND=JPRB) :: ZL_ADJ(KLON)     ! Supersat check cloud water change
REAL(KIND=JPRB) :: ZI_ADJ(KLON)     ! Supersat check cloud ice change

! required for vertical integral cloud budget
REAL(KIND=JPRB) :: ZDP(KLON)
REAL(KIND=JPRB) :: ZRHO(KLON)
REAL(KIND=JPRB) :: ZDZ(KLON)

REAL(KIND=JPRB) :: ZRTSPHY          ! Reciprocal of timestep
REAL(KIND=JPRB) :: ZSLPHY           ! Switch for SLPHY to turn off on first timestep

REAL(KIND=JPRB) :: ZEPSILON

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cloud_supersatcheck.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SLTEND',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(NDIMSLP=>YGFL%NDIMSLP, NUMFLDS=>YGFL%NUMFLDS, YA=>YGFL%YA, YL=>YGFL%YL, YI=>YGFL%YI, &
 & YO3=>YGFL%YO3, YQ=>YGFL%YQ, YR=>YGFL%YR, YS=>YGFL%YS, &
 & NCLDTOP=>YDECLDP%NCLDTOP, NSSOPT=>YDECLDP%NSSOPT, RAMIN=>YDECLDP%RAMIN, &
 & RKOOPTAU=>YDECLDP%RKOOPTAU, RTHOMO=>YDECLDP%RTHOMO, &
 & LCLDBUDC=>YDECLDP%LCLDBUDC, &
 & LCLDBUDL=>YDECLDP%LCLDBUDL, &
 & LCLDBUDI=>YDECLDP%LCLDBUDI, &
 & LCLDBUDT=>YDECLDP%LCLDBUDT, &
 & LCLDBUD_VERTINT=>YDECLDP%LCLDBUD_VERTINT, &
 & LCLDBUD_TIMEINT=>YDECLDP%LCLDBUD_TIMEINT, &
 & LBUD23=>YDEPHY%LBUD23, &
 & NSTART=>YDRIP%NSTART, RSLWX=>YDSLPHY%RSLWX, &
 & TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

ZEPSILON = 100._JPRB*EPSILON(ZEPSILON)
ZRTSPHY  = 1.0_JPRB/TSPHY

IU = 1  ! Pointer for u-wind in PSAVTEND array
IV = 2  ! Pointer for v-wind in PSAVTEND array
IT = 3  ! Pointer for temperature in PSAVTEND array

ZSLPHY = 1.0_JPRB
! If this is the first timestep then deactivate SLPHY averaging as no prev timestep
IF (NSTEP == NSTART) THEN
  ZSLPHY = 0.0_JPRB
ENDIF

!--------------------------------------------------------------------------------
!
! Semi-lagrangian averaging of physics tendencies
!
! Update prognostic variables to be averages along the semi-Lagrangian
! trajectory. Weight the "slow" tendencies to be part of this time step
! and part of the previous timestep (equally weighted if RSLWX=0.5)
!
! Note: when Yxx%LPHY false no memory is allocated for the PGFLSLP/PPHTENGFL arrays
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------- 
!
! 1. Save part of the "slow" physics tendencies (i.e. excluding the dynamics and 
!    vdf - representing the fast processes here) from this time step (t) for the next
!    time step (t+1)
!
!------------------------------------------------------------------------------- 

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

    PSAVTEND(JL,JK,IU) = PCENU(JL,JK)-PVFU(JL,JK)-PDYNU(JL,JK)
    PSAVTEND(JL,JK,IV) = PCENV(JL,JK)-PVFV(JL,JK)-PDYNV(JL,JK)
    PSAVTEND(JL,JK,IT) = PCENT(JL,JK)-PVFT(JL,JK)-PDYNT(JL,JK)-PSATT(JL,JK)
    IF (YQ%LPHY) PGFLSLP(JL,JK,YQ%MPSLP) = PCENQ(JL,JK)-PVFQ(JL,JK)-PDYNQ(JL,JK)-PSATQ(JL,JK)
    IF (YA%LPHY) PGFLSLP(JL,JK,YA%MPSLP) = PCENA(JL,JK)-PVFA(JL,JK)-PDYNA(JL,JK)-PSATA(JL,JK)
    ! Perhaps in the future we may subtract some fast tendencies (PVFCLD) as well
    IF (YL%LPHY) PGFLSLP(JL,JK,YL%MPSLP) = PCENCLD(JL,JK,NCLDQL)-PDYNCLD(JL,JK,NCLDQL)-PSATCLD(JL,JK,NCLDQL)
    IF (YI%LPHY) PGFLSLP(JL,JK,YI%MPSLP) = PCENCLD(JL,JK,NCLDQI)-PDYNCLD(JL,JK,NCLDQI)-PSATCLD(JL,JK,NCLDQI)
    IF (YR%LPHY) PGFLSLP(JL,JK,YR%MPSLP) = PCENCLD(JL,JK,NCLDQR)-PDYNCLD(JL,JK,NCLDQR)
    IF (YS%LPHY) PGFLSLP(JL,JK,YS%MPSLP) = PCENCLD(JL,JK,NCLDQS)-PDYNCLD(JL,JK,NCLDQS)
    IF (YO3%LPHY) PGFLSLP(JL,JK,YO3%MPSLP) = PCENO3(JL,JK)-PDYNO3(JL,JK)

!------------------------------------------------------------------------------- 
!
! 2. Create additional tendencies due to the averaging of this and previous time steps
!
!------------------------------------------------------------------------------- 

    PTENU(JL,JK) =  (1._JPRB-RSLWX(JK)) &
     &  *(-ZSLPHY*PSAVTEND(JL,JK,IU) + ZRTSPHY*PPHTENU(JL,JK))
    PTENV(JL,JK) =  (1._JPRB-RSLWX(JK)) &
     &  *(-ZSLPHY*PSAVTEND(JL,JK,IV) + ZRTSPHY*PPHTENV(JL,JK))
    PTENT(JL,JK) =  (1._JPRB-RSLWX(JK)) &
     &  *(-ZSLPHY*PSAVTEND(JL,JK,IT) + ZRTSPHY*PPHTENT(JL,JK))
    IF(YQ%LPHY) THEN
      PTENQ(JL,JK) =  (1._JPRB-RSLWX(JK)) &
        & *(-ZSLPHY*PGFLSLP(JL,JK,YQ%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YQ%MPSLP))
    ELSE
      PTENQ(JL,JK) = 0.0_JPRB
    ENDIF
    IF(YA%LPHY) THEN
      PTENA(JL,JK) =  (1._JPRB-RSLWX(JK)) &
        & *(-ZSLPHY*PGFLSLP(JL,JK,YA%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YA%MPSLP))
    ELSE
      PTENA(JL,JK) = 0.0_JPRB
    ENDIF
    IF(YL%LPHY) THEN
      PTENCLD(JL,JK,NCLDQL) =  (1._JPRB-RSLWX(JK)) &
       & *(-ZSLPHY*PGFLSLP(JL,JK,YL%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YL%MPSLP))
    ELSE
      PTENCLD(JL,JK,NCLDQL) = 0.0_JPRB
    ENDIF
    IF(YI%LPHY) THEN
      PTENCLD(JL,JK,NCLDQI) =  (1._JPRB-RSLWX(JK)) &
       & *(-ZSLPHY*PGFLSLP(JL,JK,YI%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YI%MPSLP))
    ELSE
      PTENCLD(JL,JK,NCLDQI) = 0.0_JPRB
    ENDIF
    IF(YR%LPHY) THEN
      PTENCLD(JL,JK,NCLDQR) =  (1._JPRB-RSLWX(JK)) &
       & *(-ZSLPHY*PGFLSLP(JL,JK,YR%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YR%MPSLP))
    ELSE
      PTENCLD(JL,JK,NCLDQR) = 0.0_JPRB
    ENDIF
    IF(YS%LPHY) THEN
      PTENCLD(JL,JK,NCLDQS) =  (1._JPRB-RSLWX(JK)) &
       & *(-ZSLPHY*PGFLSLP(JL,JK,YS%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YS%MPSLP))
    ELSE
      PTENCLD(JL,JK,NCLDQS) = 0.0_JPRB
    ENDIF
    IF(YO3%LPHY) THEN
      PTENO3(JL,JK) =  (1._JPRB-RSLWX(JK)) &
       & *(-ZSLPHY*PGFLSLP(JL,JK,YO3%MPSLP) + ZRTSPHY*PPHTENGFL(JL,JK,YO3%MPSLP))
    ELSE
      PTENO3(JL,JK) = 0.0_JPRB
    ENDIF
    
  ENDDO ! on JL
ENDDO ! on JK

!------------------------------------------------------------------------------- 
!
! 3. Perform adjustment for any supersaturation in the new state
!    Only need T,Q,A (temperature, humidity and cloud fraction)
!
!------------------------------------------------------------------------------- 

DO JK=NCLDTOP,KLEV

  DO JL=KIDIA,KFDIA
    ZT_UPD(JL) = PT(JL,JK) + TSPHY*(PCENT(JL,JK)+PTENT(JL,JK))
    ZQ_UPD(JL) = PQ(JL,JK) + TSPHY*(PCENQ(JL,JK)+PTENQ(JL,JK))
    ZA_UPD(JL) = PA(JL,JK) + TSPHY*(PCENA(JL,JK)+PTENA(JL,JK))
  ENDDO ! on JL
  
  ! Saturation adjustment step if there is any supersaturation above the defined limit
  IFTLIQICE = 1  ! Distribute liquid and ice according to mixed phase function
  CALL CLOUD_SUPERSATCHECK(YDECLDP, KIDIA, KFDIA, KLON, KLEV, IFTLIQICE, &
                         & ZT_UPD, ZQ_UPD, ZA_UPD, &
                         & PRSF1(KIDIA:KFDIA,JK), PRS1(KIDIA:KFDIA,KLEV+1), &
                         & ZT_ADJ, ZQ_ADJ, ZA_ADJ, ZL_ADJ, ZI_ADJ)

  DO JL=KIDIA,KFDIA
    ! Update the state to remove any supersaturation above defined limit
    PTENT(JL,JK) = PTENT(JL,JK) + ZT_ADJ(JL)*ZRTSPHY
    PTENQ(JL,JK) = PTENQ(JL,JK) + ZQ_ADJ(JL)*ZRTSPHY
    PTENA(JL,JK) = PTENA(JL,JK) + ZA_ADJ(JL)*ZRTSPHY
    PTENCLD(JL,JK,NCLDQL) = PTENCLD(JL,JK,NCLDQL) + ZL_ADJ(JL)*ZRTSPHY
    PTENCLD(JL,JK,NCLDQI) = PTENCLD(JL,JK,NCLDQI) + ZI_ADJ(JL)*ZRTSPHY
    ! Store for diagnostics
    ZADJA(JL,JK) = ZA_ADJ(JL)*ZRTSPHY
    ZADJL(JL,JK) = ZL_ADJ(JL)*ZRTSPHY
    ZADJI(JL,JK) = ZI_ADJ(JL)*ZRTSPHY
  ENDDO       

ENDDO ! on JK


!-----------------------------------------------------------------
! Vertical integral of all cloud process terms in one 3D field 
!-----------------------------------------------------------------

! set pointer in extra variables array
IF (LBUD23) THEN
  IS = 26  ! Take account of LBUD23 array diagnostics if turned on
ELSE
  IS = 0
ENDIF

IF (LCLDBUD_VERTINT) THEN
  DO JK=NCLDTOP,KLEV
    DO JL=KIDIA,KFDIA
 
      ZDP(JL)     = PRS1(JL,JK+1)-PRS1(JL,JK)   ! dp
      ZRHO(JL)    = PRSF1(JL,JK)/(RD*PT(JL,JK)) ! p/RT air density
      ZDZ(JL)     = ZDP(JL)/(ZRHO(JL)*RG)       ! Layer depth (m)

      PEXTRA(JL,5,1)  = PEXTRA(JL,5,1)  + ZADJA(JL,JK)*ZDZ(JL)  ! + Supersat clipping for cloud fraction
      PEXTRA(JL,20,1) = PEXTRA(JL,20,1) + ZADJL(JL,JK)*ZDZ(JL)  ! + Supersat clipping for liquid
      PEXTRA(JL,45,1) = PEXTRA(JL,45,1) + ZADJI(JL,JK)*ZDZ(JL)  ! + Supersat clipping for ice

    ENDDO
  ENDDO
  IS = IS + 1
ENDIF
!-----------------------------------------------------------------
! Cloud fraction budget 
!-----------------------------------------------------------------
IF (LCLDBUDC) THEN
  DO JK=NCLDTOP,KLEV
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZADJA(JL,JK) ! Supersat clipping
    ENDDO
  ENDDO
  IS = IS + 12
ENDIF

!-----------------------------------------------------------------
! Cloud liquid condensate budget 
!-----------------------------------------------------------------
IF (LCLDBUDL) THEN
  DO JK=NCLDTOP,KLEV
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) +  ZADJL(JL,JK) ! + Supersat clipping so far this timestep
    ENDDO
  ENDDO
  IS = IS + 22
ENDIF

!-----------------------------------------------------------------
! Cloud ice condensate budget 
!-----------------------------------------------------------------
IF (LCLDBUDI) THEN
  DO JK=NCLDTOP,KLEV
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + ZADJI(JL,JK) ! + Supersat clipping so far this timestep 
    ENDDO
  ENDDO
  IS = IS + 18
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SLTEND',1,ZHOOK_HANDLE)
END SUBROUTINE SLTEND
