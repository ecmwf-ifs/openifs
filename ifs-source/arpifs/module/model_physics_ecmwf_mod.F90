! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE MODEL_PHYSICS_ECMWF_MOD
  USE YOEPHY      , ONLY : TEPHY
  USE YOECLD      , ONLY : TECLD
  USE YOECLDP     , ONLY : TECLDP
  USE YOECND      , ONLY : TECND
  USE YOECUMF     , ONLY : TECUMF
  USE YOE_CUCONVCA, ONLY : TECUCONVCA
  USE YOEGWD      , ONLY : TEGWD
  USE YOEGWWMS    , ONLY : TEGWWMS
  USE YOETHF      , ONLY : TTHF
  IMPLICIT NONE

  TYPE MODEL_PHYSICS_ECMWF_TYPE

  TYPE(TEPHY)         :: YREPHY
  ! cloudy stuff
  TYPE(TECLD)         :: YRECLD  !! explicit clouds (diagnostic)
  TYPE(TECLDP)        :: YRECLDP !! explicit clouds (prognostic)
  TYPE(TECND)         :: YRECND  !! moist processes
  ! convection
  TYPE(TECUMF)        :: YRECUMF !! cumulus mass flux
  TYPE(TECUCONVCA)    :: YRECUCONVCA !! cellular automaton
  ! gravity wave drag
  TYPE(TEGWD)         :: YREGWD   !! parameters
  TYPE(TEGWWMS)       :: YREGWWMS !! WarnerMcIntyre GW parameterization

  TYPE(TTHF)          :: YRTHF

  CONTAINS

  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

END TYPE MODEL_PHYSICS_ECMWF_TYPE

  !---------------------------------------------------------------------

  CONTAINS

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_ECMWF_TYPE), INTENT(IN) :: SELF
  INTEGER                        , INTENT(IN) :: KDEPTH
  INTEGER                        , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_ec : '
  CALL SELF%YREPHY%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRECLD%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRECLDP%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRECND%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRECUCONVCA%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREGWD%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREGWWMS%PRINT(KDEPTH+2, KOUTNO)

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_ECMWF_MOD
