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

MODULE PTRSLB15

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     SEMI LAGRANGIAN BUFFER 15 (for adoint and tangent linear)
!     "O" stands for "origin point", "M" stands for medium point.
!     "NH" stands for "non hydrostatic", "SI" stands for "semi-implicit". 

!     * NFLDSLB15 : Number of fields in semi-lagrangian buffer 15 (TL and adjoint)
!     * RPARSL15  : parities for extrapolar extensions for interpol. buffers.
!     * MSLBUF15  : start pointer for interpolation buffer 15.

!     * MSLB1UR05 : as MSLB1UR0 but for trajectory (used in TL/AD)
!     * MSLB1VR05 : as MSLB1VR0 but for trajectory (used in TL/AD)
!     * MSLB1ZR05 : as MSLB1ZR0 but for trajectory (used in TL/AD)
!     * MSLB1WR05 : as MSLB1WR0 but for trajectory (used in TL/AD)
!     * MSLB1UR005 : as MSLB1UR00 but for trajectory (used in TL/AD)
!     * MSLB1VR005 : as MSLB1VR00 but for trajectory (used in TL/AD)
!     * MSLB1ZR005 : as MSLB1ZR00 but for trajectory (used in TL/AD)
!     * MSLB1WR005 : as MSLB1WR00 but for trajectory (used in TL/AD)
!     * MSLB1UR95 : as MSLB1UR9 but for trajectory (used in TL/AD)
!     * MSLB1VR95 : as MSLB1VR9 but for trajectory (used in TL/AD)
!     * MSLB1ZR95 : as MSLB1ZR9 but for trajectory (used in TL/AD)
!     * MSLB1U05  : as MSLB1U0 but for trajectory (used in TL/AD)
!     * MSLB1V05  : as MSLB1V0 but for trajectory (used in TL/AD)
!     * MSLB1Z05  : as MSLB1Z0 but for trajectory (used in TL/AD)
!     * MSLB1T05  : as MSLB1T0 but for trajectory (used in TL/AD)
!     * MSLB1C05  : as MSLB1C0 but for trajectory (used in TL/AD)
!     * MSLB1SP05 : as MSLB1SP0 but for trajectory (used in TL/AD)
!     * MSLB1U95  : as MSLB1U9 but for trajectory (used in TL/AD)
!     * MSLB1V95  : as MSLB1V9 but for trajectory (used in TL/AD)
!     * MSLB1Z95  : as MSLB1Z9 but for trajectory (used in TL/AD)
!     * MSLB1T95  : as MSLB1T9 but for trajectory (used in TL/AD)
!     * MSLB1GFL95: as MSLB1GFL9 but for trajectory (used in TL/AD)
!     * MSLB1C95  : as MSLB1C9 but for trajectory (used in TL/AD)
!     * MSLB1UP95 : as MSLB1UP9 but for trajectory (used in TL/AD)
!     * MSLB1VP95 : as MSLB1VP9 but for trajectory (used in TL/AD)
!     * MSLB1ZP95 : as MSLB1ZP9 but for trajectory (used in TL/AD)
!     * MSLB1TP95 : as MSLB1TP9 but for trajectory (used in TL/AD)
!     * MSLB1GFLP95: as MSLB1GFLP9 but for trajectory (used in TL/AD)
!     * MSLB1SP95 : as MSLB1SP9 but for trajectory (used in TL/AD)

TYPE :: TPTRSLB15
INTEGER(KIND=JPIM) :: NFLDSLB15
REAL(KIND=JPRB),ALLOCATABLE:: RPARSL15(:)
INTEGER(KIND=JPIM) :: MSLBUF15
INTEGER(KIND=JPIM) :: MSLB1UR05
INTEGER(KIND=JPIM) :: MSLB1VR05
INTEGER(KIND=JPIM) :: MSLB1ZR05
INTEGER(KIND=JPIM) :: MSLB1WR05
INTEGER(KIND=JPIM) :: MSLB1UR005
INTEGER(KIND=JPIM) :: MSLB1VR005
INTEGER(KIND=JPIM) :: MSLB1ZR005
INTEGER(KIND=JPIM) :: MSLB1WR005
INTEGER(KIND=JPIM) :: MSLB1UR95
INTEGER(KIND=JPIM) :: MSLB1VR95
INTEGER(KIND=JPIM) :: MSLB1ZR95
INTEGER(KIND=JPIM) :: MSLB1U05
INTEGER(KIND=JPIM) :: MSLB1V05
INTEGER(KIND=JPIM) :: MSLB1Z05
INTEGER(KIND=JPIM) :: MSLB1T05
INTEGER(KIND=JPIM) :: MSLB1C05
INTEGER(KIND=JPIM) :: MSLB1SP05
INTEGER(KIND=JPIM) :: MSLB1U95
INTEGER(KIND=JPIM) :: MSLB1V95
INTEGER(KIND=JPIM) :: MSLB1Z95
INTEGER(KIND=JPIM) :: MSLB1T95
INTEGER(KIND=JPIM) :: MSLB1GFL95
INTEGER(KIND=JPIM) :: MSLB1C95
INTEGER(KIND=JPIM) :: MSLB1UP95
INTEGER(KIND=JPIM) :: MSLB1VP95
INTEGER(KIND=JPIM) :: MSLB1ZP95
INTEGER(KIND=JPIM) :: MSLB1TP95
INTEGER(KIND=JPIM) :: MSLB1GFLP95
INTEGER(KIND=JPIM) :: MSLB1SP95
!--------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TPTRSLB15
!==============================================================

!!TYPE(TPTRSLB15), POINTER :: YRPTRSLB15 => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TPTRSLB15), INTENT(IN) :: SELF
  INTEGER         , INTENT(IN) :: KDEPTH
  INTEGER         , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%ptrslb15 : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFLDSLB15 = ', SELF%NFLDSLB15
  IF (ALLOCATED(SELF%RPARSL15)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPARSL15 allocated of shape ',SHAPE(SELF%RPARSL15),' and sum ',SUM(SELF%RPARSL15)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPARSL15 not allocated'
  ENDIF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLBUF15 = ',  SELF%MSLBUF15
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UR05 = ', SELF%MSLB1UR05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VR05 = ', SELF%MSLB1VR05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1WR05 = ', SELF%MSLB1WR05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UR95 = ', SELF%MSLB1UR95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VR95 = ', SELF%MSLB1VR95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1U05 = ', SELF%MSLB1U05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1V05 = ', SELF%MSLB1V05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1T05 = ', SELF%MSLB1T05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1C05 = ', SELF%MSLB1C05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1SP05 = ', SELF%MSLB1SP05
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1U95 = ', SELF%MSLB1U95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1V95 = ', SELF%MSLB1V95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1T95 = ', SELF%MSLB1T95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFL95 = ', SELF%MSLB1GFL95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1C95 = ', SELF%MSLB1C95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UP95 = ', SELF%MSLB1UP95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VP95 = ', SELF%MSLB1VP95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1ZP95 = ', SELF%MSLB1ZP95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1TP95 = ', SELF%MSLB1TP95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFLP95 = ', SELF%MSLB1GFLP95
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1SP95 = ', SELF%MSLB1C95

END SUBROUTINE PRINT_CONFIGURATION

END MODULE PTRSLB15
