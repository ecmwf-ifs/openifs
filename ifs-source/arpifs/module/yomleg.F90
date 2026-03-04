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

MODULE YOMLEG

USE PARKIND1 , ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Description of Legendre polynomials

! 1/ Gaussian grid and geometry:
! RW          : weights of the gaussian quadrature
! RMU         : mu              sin(theta)
! R1MU2       : 1.-MU*MU        cos(theta)**2
! R1MUI       : 1./R1MU2      1/cos(theta)**2
! R1MUA       : 1./R1MU2/RA   1/(a*cos(theta)**2)
! RSQM2       : SQRT(R1MU2)     cos(theta)
! R1QM2       : 1./SQRT(R1MU2) 1/cos(theta)
! RACTHE      : 1./SQRT(R1MU2) 1/(a*cos(theta))
! RLATIG      : arcsin(mu)      theta  GLOBAL VIEW
! RLATI       : arcsin(mu)      theta

TYPE TCSGLEG
  REAL(KIND=JPRB), ALLOCATABLE :: RW(:)
  REAL(KIND=JPRB), ALLOCATABLE :: RMU(:)
  REAL(KIND=JPRB), ALLOCATABLE :: R1MU2(:)
  REAL(KIND=JPRB), ALLOCATABLE :: R1MUI(:)
  REAL(KIND=JPRB), ALLOCATABLE :: R1MUA(:)
  REAL(KIND=JPRB), ALLOCATABLE :: RSQM2(:)
  REAL(KIND=JPRB), ALLOCATABLE :: R1QM2(:)
  REAL(KIND=JPRB), ALLOCATABLE :: RACTHE(:)
  REAL(KIND=JPRB), ALLOCATABLE :: RLATIG(:)
  REAL(KIND=JPRB), ALLOCATABLE :: RLATI(:)
END TYPE TCSGLEG


END MODULE YOMLEG
