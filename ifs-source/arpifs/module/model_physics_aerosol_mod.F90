! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 


MODULE MODEL_PHYSICS_AEROSOL_MOD
  USE YOEAERLID, ONLY : TEAERLID
  USE YOEAERMAP, ONLY : TEAERMAP
  USE YOEAERSNK, ONLY : TEAERSNK
  USE YOEAERSRC, ONLY : TEAERSRC
  USE YOEAERVOL, ONLY : TEAERVOL
  USE YOEDBUG, ONLY : TEDBUG
  IMPLICIT NONE
          
  TYPE MODEL_PHYSICS_AEROSOL_TYPE

  TYPE(TEAERLID)  :: YREAERLID  !! LIDAR simulator of 
                                !! aerosol effects
  TYPE(TEAERMAP)  :: YREAERMAP
  TYPE(TEAERSNK)  :: YREAERSNK  !! sinks
  TYPE(TEAERSRC)  :: YREAERSRC  !! sources
  TYPE(TEAERVOL)  :: YREAERVOL  !! volcanic aerosols
  TYPE(TEDBUG)    :: YREDBUG    !! aerosol debugging help

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

  END TYPE MODEL_PHYSICS_AEROSOL_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_AEROSOL_TYPE), INTENT(IN) :: SELF
  INTEGER                          , INTENT(IN) :: KDEPTH
  INTEGER                          , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_aer : '
  CALL SELF%YREAERLID%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREAERMAP%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREAERSNK%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREAERSRC%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREAERVOL%PRINT(KDEPTH+2, KOUTNO)

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_AEROSOL_MOD
