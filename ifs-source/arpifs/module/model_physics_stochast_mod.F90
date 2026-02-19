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


!!=================================

MODULE MODEL_PHYSICS_STOCHAST_MOD
  USE STOPH_MIX        , ONLY : TSTOPH
  USE YOMRANDOM_STREAMS, ONLY : TRANDOM_STREAMS
  IMPLICIT NONE

  TYPE MODEL_PHYSICS_STOCHAST_TYPE

  TYPE(TSTOPH)          :: YRSTOPH
  TYPE(TRANDOM_STREAMS) :: YR_RANDOM_STREAMS

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

  END TYPE MODEL_PHYSICS_STOCHAST_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_STOCHAST_TYPE), INTENT(IN) :: SELF
  INTEGER                           , INTENT(IN) :: KDEPTH
  INTEGER                           , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_stoch : '
  CALL SELF%YRSTOPH%PRINT(KDEPTH+2, KOUTNO)
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH+2) // 'model%yrml_phy_stoch%yr_random_streams : not yet printable'

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_STOCHAST_MOD


