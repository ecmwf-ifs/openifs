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

MODULE MODEL_ATMOS_OCEAN_COUPLING_MOD
  USE YOMMCC, ONLY : TMCC
  USE YOMCOM, ONLY : TCOM
  USE YOMCOU, ONLY : TCOU
  IMPLICIT NONE

  TYPE MODEL_ATMOS_OCEAN_COUPLING_TYPE

    TYPE(TMCC)  :: YRMCC
    TYPE(TCOM)  :: YRCOM
    TYPE(TCOU)  :: YRCOU

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
    
  END TYPE MODEL_ATMOS_OCEAN_COUPLING_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_ATMOS_OCEAN_COUPLING_TYPE), INTENT(IN) :: SELF
  INTEGER                               , INTENT(IN) :: KDEPTH
  INTEGER                               , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_aoc : '
  CALL SELF%YRMCC%PRINT(KDEPTH+2,KOUTNO)
  CALL SELF%YRCOM%PRINT(KDEPTH+2,KOUTNO)
  CALL SELF%YRCOU%PRINT(KDEPTH+2,KOUTNO)
  
  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_ATMOS_OCEAN_COUPLING_MOD
