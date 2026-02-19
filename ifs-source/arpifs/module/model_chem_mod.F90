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

MODULE MODEL_CHEM_MOD
  USE YOMOZO  , ONLY : TOZO
  USE YOMCHEM , ONLY : TCHEM
  USE YOMCOMPO, ONLY : TCOMPO
  USE DRYDEP_PAR, ONLY : TDRYDEP
  IMPLICIT NONE

  TYPE MODEL_CHEM_TYPE

    TYPE(TOZO)   :: YROZO
    TYPE(TCHEM)  :: YRCHEM
    TYPE(TCOMPO) :: YRCOMPO
    TYPE(TDRYDEP) :: YRDRYDEP

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

  END TYPE MODEL_CHEM_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_CHEM_TYPE), INTENT(IN) :: SELF
  INTEGER               , INTENT(IN) :: KDEPTH
  INTEGER               , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_chem : not yet printable'

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_CHEM_MOD
