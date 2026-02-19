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

MODULE YEMDIM

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Dimensions of model working arrays (YEMDIM)

!   NSECPLG: dimension for the Laplace operator

!*    Width of coupling/relaxation belts

!   NBZONG : half-difference between the size of C+I and C zone in meridional direction
!   NBZONL : half-difference between the size of C+I and C zone in zonal direction

!   NNOEXTZL : alternative extension zone (E') zonal dimension
!   NNOEXTZG : alternative extension zone (E') meridional dimension

!   NISNAX : zonal limit wavenumbers within the ellipse
!   NISMAX : meridional limit wavenumbers within the ellipse

!   LBIPINCI  : Boyd coupling business (...)
!   NBIPINCIX : Boyd coupling business (...)
!   NBIPINCIY : Boyd coupling business (...)

!*    Key to be moved out (in YOMOPH or so ?)

!   NEDOM  : key: options for file ALADIN :
!                        -1 --- Gridpoint file
!                         1 --- Spectral file (Aladin standard)


TYPE :: TEDIM
INTEGER(KIND=JPIM) :: NSECPLG
INTEGER(KIND=JPIM) :: NBZONG
INTEGER(KIND=JPIM) :: NBZONL
INTEGER(KIND=JPIM) :: NNOEXTZG
INTEGER(KIND=JPIM) :: NNOEXTZL
INTEGER(KIND=JPIM), POINTER :: NISMAX(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: NISNAX(:) => NULL()
LOGICAL :: LBIPINCI
INTEGER(KIND=JPIM) :: NBIPINCIX
INTEGER(KIND=JPIM) :: NBIPINCIY
INTEGER(KIND=JPIM) :: NEDOM
END TYPE TEDIM


!     ------------------------------------------------------------------
END MODULE YEMDIM
