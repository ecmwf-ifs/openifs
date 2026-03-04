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

MODULE YOMNCL

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

! ------ CLOUD CHARACTERISTICS FOR SIMPLIFIED SCHEME

! LNCLIN   : .TRUE. IF (A,L,I) grid-point upper air fields to be read
!            on input for new cloud scheme of the linearized model
! LREGCL   : .TRUE. if the regularization in the cloud scheme is used

TYPE :: TNCL
LOGICAL :: LNCLIN
LOGICAL :: LREGCL
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TNCL
!============================================================================

!!TYPE(TNCL), POINTER :: YRNCL => NULL()

!     ------------------------------------------------------------------
CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TNCL), INTENT(IN) :: SELF
INTEGER    , INTENT(IN) :: KDEPTH
INTEGER    , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_slin%yrncl : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LNCLIN = ', SELF%LNCLIN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LREGCL = ', SELF%LREGCL

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMNCL
