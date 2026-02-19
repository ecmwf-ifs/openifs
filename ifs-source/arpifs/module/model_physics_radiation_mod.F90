! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE MODEL_PHYSICS_RADIATION_MOD
  USE YOMRADF        , ONLY : TRADF
  USE YOERAD         , ONLY : TERAD
  USE YOESW          , ONLY : TESWRT
  USE YOEOVLP        , ONLY : TEOVLP
  USE YOENEUR        , ONLY : TENEUR
  USE YOELWRAD       , ONLY : TELWRAD
  USE YOEAERD        , ONLY : TEAERD
  USE YOEAERATM      , ONLY : TEAERATM
  USE YOE_UVRAD      , ONLY : TEUVRAD
  USE YOERDI         , ONLY : TERDI
  USE YOE_MCICA      , ONLY : TEMCICA
  USE YOMRCOEF       , ONLY : TRCOEF
  USE YOMTRC         , ONLY : TTRC
  USE YOMPRAD        , ONLY : RADIATION_GRID_STRUCT
  USE YOMGSGEOM      , ONLY : TGSGEOM
  USE YOERIP         , ONLY : TERIP
  USE RADIATION_SETUP, ONLY : TRADIATION
  USE EINT_MOD       , ONLY : SL_STRUCT
  USE PARKIND1       , ONLY : JPRB
  USE YOMHOOK        , ONLY : LHOOK, DR_HOOK, JPHOOK
  IMPLICIT NONE

  TYPE MODEL_PHYSICS_RADIATION_TYPE

  TYPE(TRADF)                 :: YRRADF
  TYPE(TERAD)                 :: YRERAD     !! control options for radiation config
  TYPE(TESWRT)                :: YRESWRT    !! shortwave radiation transfer coefficients
  TYPE(TEOVLP)                :: YREOVLP    !! vert distrib of cloud overlap param
  TYPE(TENEUR)                :: YRENEUR    !! neuroflux LW radiation (neuroflux radiation == coolest name ever?)
  TYPE(TELWRAD)               :: YRELWRAD   !! cloud characteristics for LW radiation
  TYPE(TEAERD)                :: YREAERD    !! spectral distribution of aerosols
  TYPE(TEAERATM)              :: YREAERATM  !! control parameters for atmos aerosols
  TYPE(TEUVRAD)               :: YREUVRAD   !! UV radiation coefs
  TYPE(TERDI)                 :: YRERDI     !! coefs in radiation interface
  TYPE(TEMCICA)               :: YREMCICA   !! cloud generator stuff
  TYPE(TRCOEF)                :: YRRCOEF    !! read & write radiation coefs
  TYPE(TTRC)                  :: YRTRC      !! storage for solar optical depths,
                                            !! perhaps only for MF phys ?
  TYPE(RADIATION_GRID_STRUCT) :: RADGRID    !! 'new'(?) radiation grid stuff
  TYPE(TGSGEOM), POINTER      :: YRAD_GSGEOM(:) => NULL()
  TYPE(TERIP)                 :: YRERIP     !! TEMPORARY TREATMENT OF RADIATION TIME, SHOULD CHANGE AT CY45
  TYPE(TRADIATION)            :: YRADIATION !! ecRad configuration
! YRRI : model grid to radiation grid
  TYPE(SL_STRUCT)             :: YRRI
! YRRO : radiation grid to model grid
  TYPE(SL_STRUCT)             :: YRRO

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

  END TYPE MODEL_PHYSICS_RADIATION_TYPE

  !---------------------------------------------------------------------

  CONTAINS

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_RADIATION_TYPE), INTENT(IN) :: SELF
  INTEGER                            , INTENT(IN) :: KDEPTH
  INTEGER                            , INTENT(IN) :: KOUTNO
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('MODEL_PHYSICS_RADIATION_MOD:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_rad : '
  CALL SELF%YRRADF%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRERAD%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRESWRT%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREOVLP%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRENEUR%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRELWRAD%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREAERD%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREAERATM%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREUVRAD%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRERDI%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREMCICA%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRRCOEF%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRTRC%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%RADGRID%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRADIATION%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRRI%PRINT(KDEPTH+2, KOUTNO, 'model%yrml_phy_rad%yrri : ')
  CALL SELF%YRRO%PRINT(KDEPTH+2, KOUTNO, 'model%yrml_phy_rad%yrro : ')
  IF (LHOOK) CALL DR_HOOK('MODEL_PHYSICS_RADIATION_MOD:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_RADIATION_MOD
