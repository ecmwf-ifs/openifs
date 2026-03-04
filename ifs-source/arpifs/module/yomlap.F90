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

MODULE YOMLAP

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Constants related to the Laplace space


!  NASN0(0:NSMAX)  : address in a spectral array of (n, m=0)
!  NASM0(0:NSMAX)  : address in a spectral array of (m, n=m)
!  NASM0G(0:NSMAX) : address in a global version of spectral array of (m, n=m)
!  NVALUE(NTPEC2)  : n value for each NSPEC2 spectral coeffient
!  MYMS(0:NUMP)    : actual wave numbers handled by this processor
!  NSPZERO(0:NSMAX): Index for the imaginary m=0 values
!  NSE0L(NUMP)     : Index array, used for example in Helmholtz operator and
!                    horizontal diffusion operator, to find the first element
!                    of the non-zero diagonals of these operators for each
!                    zonal wave number m. DM-local
!  RLAPDI(-1:NSMAX+2) :  eigen-values of the Laplace operator
!  RLAPIN(-1:NSMAX+2) :  eigen-values of its inverse

TYPE TLAP
  INTEGER(KIND=JPIM), ALLOCATABLE :: NASN0(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NASM0(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NASM0G(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NVALUE(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MYMS(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NSPZERO(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NSE0L(:)

  REAL(KIND=JPRB),    ALLOCATABLE :: RLAPDI(:)
  REAL(KIND=JPRB),    ALLOCATABLE :: RLAPIN(:)
END TYPE TLAP


END MODULE YOMLAP
