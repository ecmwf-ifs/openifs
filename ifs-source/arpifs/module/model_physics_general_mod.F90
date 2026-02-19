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

MODULE MODEL_PHYSICS_GENERAL_MOD
  USE YOMDPHY , ONLY : TDPHY
  USE YOMSLPHY, ONLY : TSLPHY
  USE YOEVDF,   ONLY : TVDF
  IMPLICIT NONE
  
  TYPE MODEL_PHYSICS_GENERAL_TYPE

    TYPE(TDPHY)  :: YRDPHY   !! dimensions
    TYPE(TSLPHY) :: YRSLPHY  !! physics & dynamics interfacing
    TYPE(TVDF)   :: YRVDF
  
  CONTAINS
  
    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

  END TYPE MODEL_PHYSICS_GENERAL_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_GENERAL_TYPE), INTENT(IN) :: SELF
  INTEGER                          , INTENT(IN) :: KDEPTH
  INTEGER                          , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_g : '
  CALL SELF%YRDPHY%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRSLPHY%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRVDF%PRINT(KDEPTH+2, KOUTNO)

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_GENERAL_MOD
