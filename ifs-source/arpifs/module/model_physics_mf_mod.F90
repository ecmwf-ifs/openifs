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

MODULE MODEL_PHYSICS_MF_MOD
  USE YOMPHY   , ONLY : TPHY
  USE YOMPHY0  , ONLY : TPHY0
  USE YOMPHY1  , ONLY : TPHY1
  USE YOMPHY2  , ONLY : TPHY2
  USE YOMPHY3  , ONLY : TPHY3
  USE YOMPHYDS , ONLY : TPHYDS
  USE YOMCVMNH , ONLY : TCVMNH
  USE YOMTOPH  , ONLY : TTOPH
  USE YOMVDOZ  , ONLY : TVDOZ
  USE YOMSIMPHL, ONLY : TSIMPHL
  USE YOMARPHY , ONLY : TARPHY
  USE YOMPARAR , ONLY : TPARAR
  USE YOMMSE   , ONLY : TMSE
  USE YOMLOUIS , ONLY : TLOUIS
  USE YOMNORGWD, ONLY : TNORGWD
  IMPLICIT NONE
  
  TYPE MODEL_PHYSICS_MF_TYPE

  TYPE(TPHY)    :: YRPHY !! (french) physics
  TYPE(TPHY0)   :: YRPHY0 !! atmospheric parameters
  TYPE(TPHY1)   :: YRPHY1 !! surface parameters
  TYPE(TPHY2)   :: YRPHY2 !! expt parameters
  TYPE(TPHY3)   :: YRPHY3 !! radiation-related parameters
  TYPE(TPHYDS)  :: YRPHYDS
  TYPE(TCVMNH)  :: YRCVMNH
  TYPE(TTOPH)   :: YRTOPH
  TYPE(TVDOZ)   :: YRVDOZ
  TYPE(TSIMPHL) :: YRSIMPHL
  TYPE(TARPHY)  :: YRARPHY
  TYPE(TPARAR)  :: YRPARAR
  TYPE(TMSE)    :: YRMSE
  TYPE(TLOUIS)  :: YRLOUIS
  TYPE(TNORGWD) :: YRNORGWD !! non-orographic GWD scheme parameters

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

  END TYPE MODEL_PHYSICS_MF_TYPE

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  USE PARKIND1  ,ONLY : JPIM
  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_MF_TYPE), INTENT(IN) :: SELF
  INTEGER(KIND=JPIM)          , INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM)          , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_mf : not yet printable'
 
  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_MF_MOD
