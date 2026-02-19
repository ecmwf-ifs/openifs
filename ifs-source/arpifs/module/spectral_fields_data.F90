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

MODULE SPECTRAL_FIELDS_DATA

!     Purpose.
!     --------
!       Type definition and data for SPECTRAL_FIELDS_MOD.
!
!     Author.
!     -------
!       Y. Tremolet
!
!     Modifications.
!     --------------
!       Original    11-Mar-2004
!       R. El Khatib 08-Jul-2014 add NPSURF in the definition of the derived type
!       T. Wilhelmsson 20-May-2015 Refactor for OOPS
! ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE
SAVE

TYPE SPECTRAL_FIELD
  REAL(KIND=JPRB), ALLOCATABLE :: SP2D(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: SP3D(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: SP1D(:,:)
  INTEGER(KIND=JPIM) :: NSMAX, NMSMAX
  INTEGER(KIND=JPIM) :: NS2D, NS3D, NS1D
  INTEGER(KIND=JPIM) :: NS2G
  INTEGER(KIND=JPIM) :: NFLEVL, NFLEVG
  INTEGER(KIND=JPIM) :: NSPEC2, NSPEC2G
  INTEGER(KIND=JPIM) :: NUMP
  INTEGER(KIND=JPIM) :: NSPSIZEL, NSPSIZEG
  INTEGER(KIND=JPIM), ALLOCATABLE :: NS2L(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NASM0(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NASM0G(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NESM0(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NESM0G(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NCPL4M(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NUMPP(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NALLMS(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MS_PROC(:,:) ! same as NALLMS, just easier to use
  INTEGER(KIND=JPIM), ALLOCATABLE :: NUMLL(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPTRLL(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPSURF(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPTRMS(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MYMS(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MYLEVS(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NGRIB(:)
  INTEGER(KIND=JPIM), POINTER :: NGRIB2(:) => NULL()
  INTEGER(KIND=JPIM), POINTER :: NGRIB3(:) => NULL()
  REAL(KIND=JPRB), POINTER :: VOR(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: DIV(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: T(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: Q(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: O3(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: L(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: I(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: SPD(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: SVD(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: NHX(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: SP(:) => NULL()
  REAL(KIND=JPRB), POINTER :: OROG(:) => NULL()
  REAL(KIND=JPRB), POINTER :: SCAL(:,:,:) => NULL() ! Scalar fields (for transforms) 
  REAL(KIND=JPRB), POINTER :: GFL(:,:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: HV(:,:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: UB(:) => NULL()
  REAL(KIND=JPRB), POINTER :: VB(:) => NULL()
END TYPE SPECTRAL_FIELD

INTEGER(KIND=JPIM) :: NPRTRV, NPRTRW, MYSETV, MYSETW
INTEGER(KIND=JPIM) :: NGRBVO, NGRBD, NGRBT, NGRBQ, NGRBO3, NGRBLNSP, NGRBCLWC, &
 & NGRBCIWC, NGRB118, NGRB119, NGRBNHX, NGRBZ
INTEGER(KIND=JPIM) :: NULOUT
LOGICAL :: LELAM

END MODULE SPECTRAL_FIELDS_DATA
