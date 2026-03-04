! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLIMAER_LAYER(YDSURF,YDMODEL,&
 &  KDIM,PAUX,STATE_T0,STATE_TMP,PSURF,GEMSL)

!**** *CLIMAER_LAYER* - Layer routine callingaerosol climatology

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities
! state_to : Derived variable for actual time model state 
! state_tmp: Derived variable for updated model state

!     ==== Input/output ====
! PSURF    : Derived variable for general surface quantities
! GEMSL    : Derived variables for local GEMS quantities


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
!      Original : 2012-12-04  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     R Hogan  3 Oct 2016: Visibility for MACC aerosol climatology
!     R Hogan 28 Mar 2017: RAD_CONFIG now in YRADIATION object
!     F Vana  14 Sep 2020: Pressure to be at t+dt
!-----------------------------------------------------------------------

USE YOECLDP                     , ONLY : NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE SURFACE_FIELDS_MIX          , ONLY : TSURF
USE TYPE_MODEL                  , ONLY : MODEL
USE PARKIND1                    , ONLY : JPIM , JPRB
USE YOMHOOK                     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER                   , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
 &                                       SURF_AND_MORE_TYPE, GEMS_LOCAL_TYPE
USE YOE_AERODIAG                , ONLY : JPAERO_WVL_AOD
USE YOMCST                      , ONLY : RG
USE RADIATION_AEROSOL_OPTICS    , ONLY : AEROSOL_SW_EXTINCTION

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                       , INTENT (INOUT) :: YDSURF
TYPE(MODEL)                       , INTENT (INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)             , INTENT (IN)    :: KDIM
TYPE (AUX_TYPE)                   , INTENT (IN)    :: PAUX
TYPE (STATE_TYPE),TARGET          , INTENT (IN)    :: STATE_T0
TYPE (STATE_TYPE),TARGET          , INTENT (IN)    :: STATE_TMP
TYPE (SURF_AND_MORE_TYPE)         , INTENT (INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)            , INTENT (INOUT) :: GEMSL
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL, JK, JAER
REAL(KIND=JPRB) :: ZDUM(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB) :: ZMACCAER(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NMCVAR) ! MACC aerosol layer mass (kg/m2)
REAL(KIND=JPRB) :: ZINV_LAYER_MASS(KDIM%KLON) ! Inverse of mass of air in a layer (m2/kg)
REAL(KIND=JPRB) :: ZMACCAER_MMR(KDIM%KLON,YDMODEL%YRML_PHY_RAD%YRADIATION%RAD_CONFIG%N_AEROSOL_TYPES) ! MACC aerosol mixing ratio (g/kg)

REAL(KIND=JPRB) :: ZTH(KDIM%KLON,KDIM%KLEV+1) ! Half level temperature (K)
REAL(KIND=JPRB) :: ZQAER(KDIM%KLON,6,KDIM%KLEV), ZTAER(KDIM%KLON,6)
REAL(KIND=JPRB) :: ZQSAT(KDIM%KLON,KDIM%KLEV) ! Saturation specific humidity (kg/kg)
REAL(KIND=JPRB) :: ZRH(KDIM%KLON)             ! Relative humidity in lowest model level

TYPE (STATE_TYPE), POINTER :: STATE
REAL(KIND=JPRB), POINTER :: ZPRES(:,:), ZPRESF(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "radaca.intfb.h"
#include "claervis.intfb.h"
#include "satur.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLIMAER_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDML_PHY_RAD=>YDMODEL%YRML_PHY_RAD,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(LAERCLIM=>YDML_PHY_RAD%YRERAD%LAERCLIM, LECSRAD=>YDML_PHY_RAD%YRERAD%LECSRAD, &
 & NCSRADF=>YDML_PHY_RAD%YRERAD%NCSRADF, YDRADIATION=>YDML_PHY_RAD%YRADIATION,&
 & YSD_VD=>YDSURF%YSD_VD, YSP_RR=>YDSURF%YSP_RR,LERADIMPL=>YDEPHY%LERADIMPL, &
 & TSPHY=>YDPHY2%TSPHY,LPHYLIN=>YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN)
!     ------------------------------------------------------------------


!*         0.     Preliminary computation

! Define state and pressure
IF(LERADIMPL) THEN
  STATE=>STATE_TMP
  ZPRES=>PAUX%PRS1
  ZPRESF=>PAUX%PRSF1
ELSE
  STATE=>STATE_T0
  ZPRES=>PAUX%PAPRS
  ZPRESF=>PAUX%PAPRSF
ENDIF


DO JK=2,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZTH(JL,JK)=(STATE%T(JL,JK-1)*ZPRESF(JL,JK-1)&
      & *(ZPRESF(JL,JK)-ZPRES(JL,JK))&
      & +STATE%T(JL,JK)*ZPRESF(JL,JK)*(ZPRES(JL,JK)-ZPRESF(JL,JK-1)))&
      & *(1.0_JPRB/(ZPRES(JL,JK)*(ZPRESF(JL,JK)-ZPRESF(JL,JK-1))))  
  ENDDO
ENDDO
DO JL=KDIM%KIDIA,KDIM%KFDIA
  ZTH(JL,1)=STATE%T(JL,1)-ZPRESF(JL,1)*(STATE%T(JL,1)-ZTH(JL,2))&
    & /(ZPRESF(JL,1)-ZPRES(JL,2))  
  ZTH(JL,KDIM%KLEV+1)=PSURF%PSP_RR(JL,YSP_RR%YT%MP9)
ENDDO

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL RADACA & CLAERVIS

! RADACA does too much for its use here: we only really need aerosol
! from one scheme at the surface yet it computes aerosols from all
! schemes at all heights. Something to improve in future.
CALL  RADACA ( YDML_PHY_RAD%YREAERD,YDML_PHY_RAD%YRERAD, YDRIP, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
   & ZPRES, PAUX%PGELAM, PAUX%PGEMU, &
   & PAUX%PCLON, PAUX%PSLON , ZTH  , &
   & ZQAER, ZMACCAER, ZDUM , GEMSL%ZCLAERS , &
   & PSURF%PSD_VD(:,YSD_VD%YAERO_WVL_DIAG(1,JPAERO_WVL_AOD)%MP), PSURF%PSD_VD(:,YSD_VD%YODSS%MP), PSURF%PSD_VD(:,YSD_VD%YODDU%MP),&
   & PSURF%PSD_VD(:,YSD_VD%YODOM%MP), PSURF%PSD_VD(:,YSD_VD%YODBC%MP), PSURF%PSD_VD(:,YSD_VD%YODSU%MP))

IF ((YDML_PHY_RAD%YREAERATM%LAERCCN .OR. YDML_PHY_RAD%YREAERATM%LAERRRTM .OR. YDML_PHY_RAD%YRERAD%NAERMACC == 1) &
     &  .AND. .NOT. YDML_PHY_RAD%YRERAD%LUSEPRE2017RAD) THEN
  ! We are using the MACC/CAMS aerosol climatology AND the newer
  ! radiation scheme: overwrite the extinction coefficient in
  ! GEMSL%ZCLAERS with the value for the MACC/CAMS aerosol climatology
  CALL SATUR(KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, KDIM%KLEV, LPHYLIN, &
       &  ZPRESF, STATE%T, ZQSAT, 2)
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZRH(JL) = STATE%Q(JL,KDIM%KLEV)/ZQSAT(JL,KDIM%KLEV)
    ZINV_LAYER_MASS(JL) = RG / (ZPRESF(JL,KDIM%KLEV)-ZPRESF(JL,KDIM%KLEV-1))
  ENDDO

  ! Convert aerosol layer mass in kg/m2 to mass mixing ratio in g/kg
  !
  ! NOTE: This whole mechanism appears to be predicated on the assumption that the set of
  !   aerosol types used here (unconditionally from the climatology?) is the same as that
  !   used in the radiation scheme (from prognostics if LAERRRTM=true).
  !
  !   If radiation is using prognostic aerosol but visibility is using the climatology,
  !   don't we then really need separate instances of RAD_CONFIG?
  !
  !   I've patched it up to not actually crash at least when N_AEROSOL_TYPES /= NMCVAR, but
  !   this will only behave correctly if the types in the climatology match the first
  !   NMCVAR types in the prognostic configuration...
  !
  !   Zak 2018-11-12
  DO JAER = 1,MIN(YDML_PHY_RAD%YRERAD%NMCVAR,YDRADIATION%RAD_CONFIG%N_AEROSOL_TYPES)
    DO JL = KDIM%KIDIA,KDIM%KFDIA
      ZMACCAER_MMR(JL,JAER) = 1000.0_JPRB*ZMACCAER(JL,KDIM%KLEV,JAER)*ZINV_LAYER_MASS(JL)
    ENDDO
  ENDDO
  IF (YDRADIATION%RAD_CONFIG%N_AEROSOL_TYPES > YDML_PHY_RAD%YRERAD%NMCVAR) THEN
    ZMACCAER_MMR(KDIM%KIDIA:KDIM%KFDIA,YDML_PHY_RAD%YRERAD%NMCVAR+1:YDRADIATION%RAD_CONFIG%N_AEROSOL_TYPES) = 0._JPRB
  ENDIF


  ! Use new radiation scheme to compute extinction coefficient in the
  ! 10th shortwave band, which includes 550 nm; by providing mass
  ! mixing ratio in g/kg we get output in km-1, which is what CLAERVIS
  ! expects. This non-SI stuff wastes everyone's time and should be
  ! fixed!
  CALL AEROSOL_SW_EXTINCTION(KDIM%KLON, KDIM%KIDIA, KDIM%KFDIA, &
       &  YDRADIATION%RAD_CONFIG, 10, ZMACCAER_MMR, ZRH, GEMSL%ZCLAERS)

ENDIF

! Calculate horizontal visibility 
CALL  CLAERVIS ( KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
     & GEMSL%ZCLAERS, STATE%A, STATE%CLD(:,:,NCLDQL), STATE%CLD(:,:,NCLDQI),&
     & STATE%CLD(:,:,NCLDQR),  STATE%CLD(:,:,NCLDQS), ZPRESF, STATE%T, &
     & GEMSL%ZVISICL )     

!-- for diagnostics only: store the 2D distribution of climatological aerosols in EXTRA fields    
IF (LECSRAD .AND. NCSRADF == 1 .AND. LAERCLIM) THEN
  DO JAER=1,6
    ZTAER(KDIM%KIDIA:KDIM%KFDIA,JAER)=0._JPRB        
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA    
        ZTAER(JL,JAER)=ZTAER(JL,JAER)+ZQAER(JL,JAER,JK)
      ENDDO
    ENDDO
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PSURF%PSD_XA(JL,JAER+8,14)=PSURF%PSD_XA(JL,JAER+8,14)+ZTAER(JL,JAER)*TSPHY
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLIMAER_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLIMAER_LAYER
