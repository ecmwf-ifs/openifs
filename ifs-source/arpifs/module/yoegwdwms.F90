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

MODULE YOEGWDWMS

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

! ------  CONTROLS FOR SIMPLIFIED NON-OROGRAPHIC GRAVITY WAVE SCHEME

! LREGNOGWD: .TRUE. if the regularization for non-orographic GWD is used

TYPE :: TEGWDWMS
LOGICAL :: LREGWWMS
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEGWDWMS
!============================================================================

!!TYPE(TEGWDWMS), POINTER :: YREGWDWMS => NULL()

!     --------------------------------------------------------------------
CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TEGWDWMS), INTENT(IN) :: SELF
  INTEGER        , INTENT(IN) :: KDEPTH
  INTEGER        , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH)    // 'model%yrml_phy_slin%yregwdwms : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LREGWWMS = ', SELF%LREGWWMS

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEGWDWMS
