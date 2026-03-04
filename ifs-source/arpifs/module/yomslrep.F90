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

MODULE YOMSLREP

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     USED FOR BIT REPRODUCIBILITY OF SEMI-LAGRANGIAN ADJOINT (SLAD)

!     NADMAP     : USED FOR REMAPPING POINTS IN SL BUFFER SO THAT THEY
!                : ARE UNIQUE
!     NADCORE    : INNER HALO FOR ADJOINT INTERPOLATIONS
!                : ARE UNIQUE
!     LADCORE    : T - IF ONE OF THE NGPTOTAD CORE POINTS, F OTHERWISE
!     NGPTOTAD   : NUMBER OF GRID POINTS IN INNER HALO
!     RSASIGN    : USED FOR CORRECTING SYMMETRIC / ANTISYMMETRIC PROPERTIES
!                : (MIRRORED LATITUDES NEAR POLES)

TYPE :: TSLREP
INTEGER(KIND=JPIM) :: NGPTOTAD

INTEGER(KIND=JPIM),ALLOCATABLE :: NADMAP(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NADCORE(:)

LOGICAL,ALLOCATABLE :: LADCORE(:)

REAL(KIND=JPRB),ALLOCATABLE :: RSASIGN(:,:)
 
!     ------------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TSLREP
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------


!!TYPE(TSLREP), POINTER :: YRSLREP => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TSLREP), INTENT(IN) :: SELF
  INTEGER      , INTENT(IN) :: KDEPTH
  INTEGER      , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_dyn%yrslrep : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGPTOTAD = ', SELF%NGPTOTAD
  IF (ALLOCATED(SELF%NADMAP)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NADMAP allocated of shape ', SHAPE(SELF%NADMAP)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NADMAP not allocated'
  ENDIF
  IF (ALLOCATED(SELF%NADCORE)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NADCORE allocated of shape ', SHAPE(SELF%NADCORE)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NADCORE not allocated'
  ENDIF
  IF (ALLOCATED(SELF%LADCORE)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LADCORE allocated of shape ', SHAPE(SELF%LADCORE),' values : ',SELF%LADCORE
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LADCORE not allocated'
  ENDIF

  IF (ALLOCATED(SELF%RSASIGN)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSASIGN allocated of shape ', SHAPE(SELF%RSASIGN),' sum : ',SUM(SELF%RSASIGN)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSASIGN not allocated'
  ENDIF
  
END SUBROUTINE PRINT_CONFIGURATION
!     ------------------------------------------------------------------
END MODULE YOMSLREP
