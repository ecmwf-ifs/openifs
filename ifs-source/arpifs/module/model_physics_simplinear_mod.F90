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

MODULE MODEL_PHYSICS_SIMPLINEAR_MOD
  USE YOEPHLI    , ONLY : TEPHLI
  USE YOMCUMFS   , ONLY : TCUMFS
  USE YOEGWDWMS  , ONLY : TEGWDWMS
  USE YOECUMF2   , ONLY : TECUMF2
  USE YOPHLC     , ONLY : TPHLC
  USE YOPHNC     , ONLY : TPHNC
  USE YOMNCL     , ONLY : TNCL
  USE YOMSRFTLAD , ONLY : TSRFTLAD
  USE YOMSPHYHIST, ONLY : SPHYS_HIST_TYPE
  IMPLICIT NONE

  TYPE MODEL_PHYSICS_SIMPLINEAR_TYPE

  TYPE(TEPHLI)   :: YREPHLI   !! constants for lin-phys
  TYPE(TCUMFS)   :: YRCUMFS   !! simplified convection scheme
  TYPE(TEGWDWMS) :: YREGWDWMS !! simplified non-oro g-w scheme
  TYPE(TECUMF2)  :: YRECUMF2  !! new linearized cumulus mass flux
  TYPE(TPHLC)    :: YRPHLC    !! simplified phys switches
  TYPE(TPHNC)    :: YRPHNC    !! switches for t-dt trajectory & phys
  TYPE(TNCL)     :: YRNCL     !! simplified cloud characteristics
  TYPE(TSRFTLAD) :: YRSRFTLAD !! lin phys skin temp
  TYPE(SPHYS_HIST_TYPE), POINTER :: GPHIST(:) => NULL()  !! trajectory “memory”
  !---------------------------------------------------------------------
  
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
                                                       !! across timesteps
  END TYPE MODEL_PHYSICS_SIMPLINEAR_TYPE
  !======================================================================

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(MODEL_PHYSICS_SIMPLINEAR_TYPE), INTENT(IN) :: SELF
  INTEGER                             , INTENT(IN) :: KDEPTH
  INTEGER                             , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_phy_slin : '
  CALL SELF%YREPHLI%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRCUMFS%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YREGWDWMS%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRECUMF2%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRPHLC%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRPHNC%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRNCL%PRINT(KDEPTH+2, KOUTNO)
  CALL SELF%YRSRFTLAD%PRINT(KDEPTH+2, KOUTNO)
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH+2) // 'model%yrml_phy_slin%gphist not printed as it is trajectory storage'

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_PHYSICS_SIMPLINEAR_MOD
