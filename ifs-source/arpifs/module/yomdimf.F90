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

MODULE YOMDIMF

USE PARKIND1 , ONLY : JPIM

IMPLICIT NONE

SAVE

TYPE :: TDIMF

! === NUMBER OF FIELDS ========================================================

! NFTHER : number of spectral thermodynamic variables
! NF3D   = number of 3D fields in the state of the model
! NFD2D  : number of 2D fields in the dynamics
! NFC2D  : number of 2D fields in the boundaries
! NPPM   : Number of interpolation methods in post-processing
! NS3D   : number of 3D fields in spectral space
! NS2D   : number of 2D fields in spectral space
! NS1D   : number of 1D fields in spectral space (for Aladin consistency -CF)
! NGRBSP2, NGRBSP3 - Grib codes for fields in SPA2 and SPA3

INTEGER(KIND=JPIM) :: NFTHER
INTEGER(KIND=JPIM) :: NF3D
INTEGER(KIND=JPIM) :: NFD2D
INTEGER(KIND=JPIM) :: NFC2D
!!INTEGER(KIND=JPIM) :: NPPM  !! moved to PARDIM in CY45
INTEGER(KIND=JPIM) :: NS3D
INTEGER(KIND=JPIM) :: NS2D
INTEGER(KIND=JPIM) :: NS1D
INTEGER(KIND=JPIM),ALLOCATABLE :: NGRBSP3(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NGRBSP2(:)


! === GRID POINT ARRAYS =======================================================

! LVOR  : controls the allocation of vorticity
! LADER : controls the allocation of vor div and derivatives
! LUVDER: controls the allocation of derivatives for u and v
! LSPT  : .TRUE. if temperature variable as spectral field

LOGICAL :: LVOR
LOGICAL :: LADER
LOGICAL :: LUVDER
LOGICAL :: LSPT
 
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TDIMF

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TDIMF), INTENT(IN) :: SELF
  INTEGER     , INTENT(IN) :: KDEPTH
  INTEGER     , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%yrdimf : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFTHER = ', SELF%NFTHER
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NF3D = ', SELF%NF3D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFD2D = ', SELF%NFD2D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFC2D = ', SELF%NFC2D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NS3D = ', SELF%NS3D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NS2D = ', SELF%NS2D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NS1D = ', SELF%NS1D
  IF (ALLOCATED(SELF%NGRBSP3)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGRBSP3 allocated, size: ', SIZE(SELF%NGRBSP3)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGRBSP3 not allocated'
  END IF
  IF (ALLOCATED(SELF%NGRBSP2)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGRBSP2 allocated, size: ', SIZE(SELF%NGRBSP2)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGRBSP2 not allocated'
  END IF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVOR = ', SELF%LVOR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LADER = ', SELF%LADER
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LUVDER = ', SELF%LUVDER
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPT = ', SELF%LSPT
  WRITE(KOUTNO,*) ''

END SUBROUTINE PRINT_CONFIGURATION


!!TYPE(TDIMF), POINTER :: YRDIMF => NULL()

!     ------------------------------------------------------------------
END MODULE YOMDIMF
