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

MODULE YOMTNH

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!   -----------------------------------------------------------------

!*   Buffer for trajectory array at time t-dt (nonhydrostatic dynamics)

!    NLENNH95B                    : length of buffer for 3D fields
!    TRAJNH(NLENNH95B,0:NSTOP)    : buffer for 3D fields

TYPE :: TTNH
INTEGER(KIND=JPIM) :: NLENNH95B
REAL(KIND=JPRB),ALLOCATABLE:: TRAJNH(:,:)
!--------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TTNH

!!TYPE(TTNH), POINTER :: YRTNH => NULL()
!==============================================================

!!TYPE(TPTRSLB1), POINTER :: YRPTRSLB1 => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TTNH), INTENT(IN) :: SELF
  INTEGER    , INTENT(IN) :: KDEPTH
  INTEGER    , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%yrtnh : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLENNH95B = ', SELF%NLENNH95B
  IF (ALLOCATED(SELF%TRAJNH)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRAJNH allocated of shape ', SHAPE(SELF%TRAJNH),' sum ',SUM(SELF%TRAJNH)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRAJNH not allocated'
  ENDIF

END SUBROUTINE PRINT_CONFIGURATION
!   ----------------------------------------------------------------
END MODULE YOMTNH
