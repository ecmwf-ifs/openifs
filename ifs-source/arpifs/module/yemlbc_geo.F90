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

MODULE YEMLBC_GEO

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Defining specific geometry variables for LAM model: those linked with LBC and I+E zones.
!     These variables are set-up in SUEGEOLBC and should not be modified elsewhere.

!     ------------------------------------------------------------------

TYPE :: TELBC_GEO

  ! NEDLST         : Nb. of points in E+I area (DM-local var.)
  ! NIND_LIST,NIND_LEN: help arrays for memory transfers between
  !  (NPROMA,NGPBLKS)-dimensioned arrays and NEDLST-dimensioned arrays.

  INTEGER(KIND=JPIM) :: NEDLST
  INTEGER(KIND=JPIM), ALLOCATABLE :: NIND_LIST(:,:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NIND_LEN(:)

END TYPE TELBC_GEO

!TYPE(TELBC_GEO), POINTER :: YRGEOLBC => NULL()   ! moved to type_geometry.F90

!     ------------------------------------------------------------------
END MODULE YEMLBC_GEO
