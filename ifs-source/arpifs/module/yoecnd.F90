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

MODULE YOECND

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------

TYPE :: TECND
REAL(KIND=JPRB) :: REPFLM
REAL(KIND=JPRB) :: REPQMI
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TECND
!============================================================================

!!TYPE(TECND), POINTER :: YRECND => NULL()

!     -----------------------------------------------------------------
!*    CONTROL PARAMETERS FOR MOIST PROCESSES

! REPFLM :  Minimum flux to avoid zero division in ice proportion
!           computations
! REPQMI :  Minimum specific humidity (security within QNEGAT)
!     -----------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TECND), INTENT(IN) :: SELF
  INTEGER     , INTENT(IN) :: KDEPTH
  INTEGER     , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_ec%yrecnd : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPFLM = ',SELF%REPFLM
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPQMI = ',SELF%REPQMI
 
END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOECND
