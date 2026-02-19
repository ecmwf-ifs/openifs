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

MODULE YOEOVLP

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOEOVLP* VERTICAL DISTRIBUTION OF CLOUD OVERLAP PARAMETER
!     ------------------------------------------------------------------

TYPE :: TEOVLP
REAL(KIND=JPRB),ALLOCATABLE:: RA1OVLP(:)
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEOVLP
!============================================================================

!!TYPE(TEOVLP), POINTER :: YREOVLP => NULL()

!     J.-J. MORCRETTE    E.C.M.W.F.     01/02/16

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *RA1OVLP* REAL      Alpha1 (Hogan, Illingworth, 2001)

!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEOVLP), INTENT(IN) :: SELF
INTEGER      , INTENT(IN) :: KDEPTH
INTEGER      , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yreovlp : '
IF (ALLOCATED(SELF%RA1OVLP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RA1OVLP ALLOCATED OF SHAPE ', SHAPE(SELF%RA1OVLP), &
 &        ' SUM ',SUM(SELF%RA1OVLP)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEOVLP
