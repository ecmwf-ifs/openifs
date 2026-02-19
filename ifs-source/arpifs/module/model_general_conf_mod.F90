! (C) Copyright 2017- ECMWF.
! (C) Copyright 2017- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE MODEL_GENERAL_CONF_MOD
  USE TYPE_GEOMETRY, ONLY : GEOMETRY
  USE YOMDIMF      , ONLY : TDIMF
  USE YOM_YGFL     , ONLY : TYPE_GFLD
  USE YOMRIP       , ONLY : TRIP
  USE YOMMODERRMOD , ONLY : TMODERR
  USE YOMSPSDT     , ONLY : TSPPT_CONFIG
  USE SPP_MOD      , ONLY : TSPP_CONFIG
  USE TYPE_ECV     , ONLY : TECVDIM
  USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
  IMPLICIT NONE

  TYPE MODEL_GENERAL_CONF_TYPE

    TYPE(GEOMETRY), POINTER :: GEOM => NULL()

    TYPE(TDIMF)             :: YRDIMF                  !! number of fields
    TYPE(TYPE_GFLD)         :: YGFL                    !! gfl descriptors
    TYPE(TRIP)              :: YRRIP                   !! TEMPORARY TREATMENT OF TIME, SHOULD CHANGE AT CY45
    TYPE(TECVDIM)           :: YRDIMECV                !! ECV field
    TYPE(TMODERR)           :: YRMODERR                !! Model error config
    TYPE(TSPPT_CONFIG)      :: YRSPPT_CONFIG           !! SPPT config
    TYPE(TSPP_CONFIG)       :: YRSPP_CONFIG            !! SPP config
    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

  END TYPE MODEL_GENERAL_CONF_TYPE

  !---------------------------------------------------------------------

  CONTAINS

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

    IMPLICIT NONE
    CLASS(MODEL_GENERAL_CONF_TYPE), INTENT(IN) :: SELF
    INTEGER                       , INTENT(IN) :: KDEPTH
    INTEGER                       , INTENT(IN) :: KOUTNO

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODEL_GENERAL_CONF_MOD:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)

    WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_gconf : '
    CALL SELF%YRDIMF%PRINT(KDEPTH+2,KOUTNO)
    CALL SELF%YRRIP%PRINT(KDEPTH+2,KOUTNO)
    CALL SELF%YRMODERR%PRINT(KDEPTH+2,KOUTNO)

    IF (LHOOK) CALL DR_HOOK('MODEL_GENERAL_CONF_MOD:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_GENERAL_CONF_MOD
