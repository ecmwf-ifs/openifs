! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADVIS_LAYER(YDMODEL,YDERAD,YGFL,KDIM,PAUX,STATE,GEMSL,PDIAG)

!**** *RADVIS_LAYER* - Layer routine calling visibility diagnostics

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities
! STATE    : Derived variable for model state

!     ==== Input/output ====
! GEMSL    : Derived variable for local GEMS quantities
! PDIAG    : Derived variable for diagn. quantities


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

!     AUTHOR.
!     -------
!      Original : 11-Feb-2012  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      R. Hogan   25-Nov-2022  Adapt to mixtures of clim./prog. aerosol & replace functionality in CLIMAER_LAYER
  
!-----------------------------------------------------------------------

USE TYPE_MODEL,ONLY : MODEL
USE YOECLDP  , ONLY : NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE YOERAD   , ONLY : TERAD
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &                   GEMS_LOCAL_TYPE, AUX_DIAG_TYPE
USE YOM_YGFL , ONLY : TYPE_GFLD
USE RADIATION_AEROSOL_OPTICS    , ONLY : AEROSOL_EXTINCTION
USE YOEAERRADDESC               , ONLY : IPROGNOSTIC

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (MODEL)                   , INTENT(INOUT) :: YDMODEL
TYPE (TERAD)                   , INTENT(INOUT) :: YDERAD
TYPE (TYPE_GFLD)               , INTENT(INOUT) :: YGFL
TYPE (DIMENSION_TYPE)          , INTENT(IN)    :: KDIM
TYPE (AUX_TYPE)                , INTENT(IN)    :: PAUX
TYPE (STATE_TYPE)              , INTENT(IN)    :: STATE
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL, JK

 ! Local value of aerosol extinction coefficient (m-1)
REAL(KIND=JPRB) :: ZAEREXTSURF(KDIM%KLON)

! Aerosol mixing ratio in lowest layer(kg/kg)
REAL(KIND=JPRB) :: ZAER_MMR(KDIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%YAER_RAD_DESC%NAEROSOL)

REAL(KIND=JPRB) :: ZRH(KDIM%KLON)             ! Relative humidity in lowest model level
REAL(KIND=JPRB) :: ZQSAT(KDIM%KLON,KDIM%KLEV) ! Saturation specific humidity (kg/kg)

! Used in call to RADACA: if we remove RADACA, remove these
REAL(KIND=JPRB) :: ZTH(KDIM%KLON,KDIM%KLEV+1) ! Half level temperature (K)
REAL(KIND=JPRB) :: ZQAER(KDIM%KLON,6,KDIM%KLEV)
REAL(KIND=JPRB) :: ZDUM(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB) :: ZMACCAER(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%YAER_RAD_DESC%NAEROSOL) ! aerosol layer mass (kg/m2)
REAL(KIND=JPRB) :: ZODTO(KDIM%KLON), ZODSS(KDIM%KLON), ZODDU(KDIM%KLON) ! Dummy outputs...
REAL(KIND=JPRB) :: ZODOM(KDIM%KLON), ZODBC(KDIM%KLON), ZODSU(KDIM%KLON) ! ...from RADACA

! Pressure on half and full levels (Pa)
REAL(KIND=JPRB), POINTER :: ZPRES(:,:), ZPRESF(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "radvis.intfb.h"
#include "satur.intfb.h"
#include "merge_aer_surface.intfb.h"
#include "radaca.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RADVIS_LAYER',0,ZHOOK_HANDLE)

ASSOCIATE(NACTAERO    => YGFL%NACTAERO, &
     &    YREAERATM   => YDMODEL%YRML_PHY_RAD%YREAERATM, &
     &    YDRADIATION => YDMODEL%YRML_PHY_RAD%YRADIATION, &
     &    YDML_PHY_RAD=> YDMODEL%YRML_PHY_RAD, &
     &    YDRIP       => YDMODEL%YRML_GCONF%YRRIP)
  
!     ------------------------------------------------------------------

!*         0.     Preliminaries

! Define state and pressure
IF(YDMODEL%YRML_PHY_EC%YREPHY%LERADIMPL) THEN
  ! Next timestep
  ZPRES  => PAUX%PRS1
  ZPRESF => PAUX%PRSF1
ELSE
  ! Current timestep
  ZPRES  => PAUX%PAPRS
  ZPRESF => PAUX%PAPRSF
ENDIF


!*         1.     Compute or copy aerosol extinction coefficient in lowest model layer

IF (NACTAERO > 0 .AND. YREAERATM%LPROGAERVIS) THEN
  ! We have prognostic aerosols and are using them in visibility
  ! calculation: copy over extinction coefficient

  DO JL = KDIM%KIDIA,KDIM%KFDIA
    ZAEREXTSURF(JL) = GEMSL%AEREXTSURF(JL)
  ENDDO
  
ELSEIF ((YREAERATM%LAERCCN .OR. YREAERATM%LAERRRTM .OR. YDERAD%NAERMACC == 1) &
     &  .AND. .NOT. YDERAD%LUSEPRE2017RAD) THEN
  ! Use the aerosols seen by the ecRad radiation scheme, which may be
  ! prognostic, climatological or a mixture

  ZAER_MMR = 0.0_JPRB

  ! Add surface climatological aerosols to ZAER_MMR as mass mixing
  ! ratios in kg/kg
  IF (YDERAD%YAER_RAD_DESC%NCLIMATOLOGICAL > 0) THEN
    CALL YDRIP%YREAERC%CALC_SURFACE(KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV, &
         &                    PAUX%PGEMU,PAUX%PGELAM,ZPRES,ZAER_MMR,LDLEAVESPACE=.TRUE.)
  ENDIF

  ! Merge in the prognostic aerosols
  IF (YDERAD%YAER_RAD_DESC%NPROGNOSTIC > 0) THEN
    CALL MERGE_AER_SURFACE(YDERAD%YAER_RAD_DESC,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON, &
         &  IPROGNOSTIC, GEMSL%ZAEROP(:,KDIM%KLEV,:), ZAER_MMR, LDCONTIG=.FALSE.)
  ENDIF

  ! Compute relative humidity
  CALL SATUR(YDMODEL%YRML_PHY_EC%YRTHF, YDMODEL%YRCST, KDIM%KIDIA, KDIM%KFDIA, &
    & KDIM%KLON, KDIM%KLEV, KDIM%KLEV, YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
    &  ZPRESF, STATE%T, ZQSAT, 2)
  DO JL = KDIM%KIDIA,KDIM%KFDIA
    ZRH(JL) = STATE%Q(JL,KDIM%KLEV)/ZQSAT(JL,KDIM%KLEV)
  ENDDO

  ! Use ecRad to compute extinction coefficient at 550 nm (m-1)
  CALL AEROSOL_EXTINCTION(KDIM%KLON, KDIM%KIDIA, KDIM%KFDIA, &
       &  YDRADIATION%RAD_CONFIG, 550.0e-9_JPRB, ZAER_MMR, ZRH, ZAEREXTSURF)
  
ELSE
  ! Use the Tegen climatology: first calculate half-level temperatures
  DO JK = 2,KDIM%KLEV
    DO JL = KDIM%KIDIA,KDIM%KFDIA
      ! Pressure-weighted linear interpolation in pressure
      ZTH(JL,JK) = (STATE%T(JL,JK-1)*ZPRESF(JL,JK-1)&
           & *(ZPRESF(JL,JK)-ZPRES(JL,JK))&
           & +STATE%T(JL,JK)*ZPRESF(JL,JK)*(ZPRES(JL,JK)-ZPRESF(JL,JK-1)))&
           & *(1.0_JPRB/(ZPRES(JL,JK)*(ZPRESF(JL,JK)-ZPRESF(JL,JK-1))))  
    ENDDO
  ENDDO
  DO JL = KDIM%KIDIA,KDIM%KFDIA
    ! Clamp top temperature
    ZTH(JL,1) = STATE%T(JL,1)
    ! Extrapolate bottom temperature
    ZTH(JL,KDIM%KLEV+1) = 2.0_JPRB * STATE%T(JL,KDIM%KLEV) &
         &  - ZTH(JL,KDIM%KLEV)
  ENDDO

  ! Use Tegen aerosols, returned in
  ! ZAEREXTSURF - all the other outputs are ignored. Note that this
  ! call is highly inefficient: RADACA does much more than compute
  ! surface extinction coefficient.
  CALL  RADACA ( YDML_PHY_RAD%YREAERD,YDML_PHY_RAD%YRERAD, YDRIP, &
       & KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
       & ZPRES, PAUX%PGELAM, PAUX%PGEMU, &
       & PAUX%PCLON, PAUX%PSLON , ZTH  , &
       & ZQAER, ZMACCAER, ZDUM , ZAEREXTSURF, &
       & ZODTO, ZODSS, ZODDU, ZODOM, ZODBC, ZODSU)

ENDIF

!*         2.     Compute visibility

CALL RADVIS(KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV , &
     &  PAUX%PRSF1, STATE%T, STATE%A, STATE%CLD(:,:,NCLDQL), STATE%CLD(:,:,NCLDQI), &
     &  STATE%CLD(:,:,NCLDQR), STATE%CLD(:,:,NCLDQS), &
     &  ZAEREXTSURF, PDIAG%PVISIH)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADVIS_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE RADVIS_LAYER
