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

MODULE YOMCOU

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

!*    INFORMATION FOR USE WITH OASIS COUPLER

! NPIAT      PROCESS IDENTIFIER OF THE ATMOSPHERE RUN
! NCULMR     C UNIT IDENTIFIER FOR OASIS MODEL READ PIPE
! NCULMW     C UNIT IDENTIFIER FOR OASIS MODEL WRITE PIPE
! NCULF(20)  C UNIT IDENTIFIERS FOR OASIS MODEL FIELD PIPES

TYPE :: TCOU
INTEGER(KIND=JPIM) :: NCULF(0:20)
INTEGER(KIND=JPIM) :: NPIAT
INTEGER(KIND=JPIM) :: NCULMR
INTEGER(KIND=JPIM) :: NCULMW
 !---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE  TCOU
!======================================================================

!!TYPE(TCOU), POINTER :: YRCOU => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TCOU), INTENT(IN) :: SELF
INTEGER    , INTENT(IN) :: KDEPTH
INTEGER    , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC 

IDEPTHLOC = KDEPTH + 2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_aoc%yrcou : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCULF sum = ', SUM(SELF%NCULF)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPIAT = ', SELF%NPIAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCULMR = ', SELF%NCULMR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCULMW = ', SELF%NCULMW

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMCOU
