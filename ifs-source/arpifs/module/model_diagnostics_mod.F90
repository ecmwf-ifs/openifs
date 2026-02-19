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

MODULE MODEL_DIAGNOSTICS_MOD
  USE YOMCDDH , ONLY : TCDDH
  USE YOMLDDH , ONLY : TLDDH
  USE YOMMDDH , ONLY : TMDDH
  USE YOMSDDH , ONLY : TSDDH
  USE YOMTDDH , ONLY : TTDDH
  USE YOMGPDDH, ONLY : TGPDDH
  USE YOMPADDH, ONLY : TPADDH
  USE YOMSPDDH, ONLY : TSPDDH
  IMPLICIT NONE

  TYPE MODEL_DIAGNOSTICS_TYPE

!!  TYPE(GRIDPOINT_BUFFER)   :: GFUBUF, XFUBUF, GPPCBUF
  TYPE(TCDDH)              :: YRCDDH
  TYPE(TLDDH)              :: YRLDDH
  TYPE(TMDDH)              :: YRMDDH
  TYPE(TSDDH)              :: YRSDDH
  TYPE(TTDDH)              :: YRTDDH
  TYPE(TGPDDH)             :: YRGPDDH
  TYPE(TPADDH)             :: YRPADDH
  TYPE(TSPDDH)             :: YRSPDDH

CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
    
  END TYPE MODEL_DIAGNOSTICS_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_DIAGNOSTICS_TYPE), INTENT(IN) :: SELF
  INTEGER                      , INTENT(IN) :: KDEPTH
  INTEGER                      , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_diag : not yet printable'

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_DIAGNOSTICS_MOD
