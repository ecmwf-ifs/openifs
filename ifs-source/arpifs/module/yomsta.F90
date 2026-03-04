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

MODULE YOMSTA

USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Description of standard atmosphere (prefix R)
!*    and of reference atm. used in post-processing (prefix V)

!    Suffixe    for RZ, RT and RP (VZ, VT AND VP)            |      RDTDZ
!                                                            |
!  1:(SUR)   Soil Surface                                    | 1: Troposphere
!  2:(TROP)  Bottom of tropopause                            | 2: Tropopause
!  3:(STRA)  Top of tropopause - Bottom of stratosphere      | 3: Stratosphere
!  4:(STR2)  Modification of temperature gradient            | 4:      "
!  5:(STPO)  Top of stratosphere - Bottom of stratopause     | 5: Stratopause
!  6:(MESO)  Top of stratopause - Bottom of mesosphere       | 6: Mesosphere
!  7:(MES2)  Modification of temperature gradient            | 7:      "
!  8:(MEPO)  Top of mesosphere - Bottom of mesopause         | 8: Mesopause
!  9:(ABOV)  Top of mesopause                                | 9: Above

!  and we have for example :
!     RDTDZ1: DT/DZ below tropopause (-6.5 K/KM)
!     RZTROP: Top of troposphere (Height of tropopause : 11000 M)
!     RTSUR : Surface temperature (288.15)

!     STPREH(0:NFLEVG) : PRESSURE OF THE HALF LEVELS
!     STPRE   (NFLEVG) : PRESSURE OF THE LEVELS
!     STPHI   (NFLEVG) : GEOPOTENTIAL
!     STTEM   (NFLEVG) : TEMPERATURE
!     STDEN   (NFLEVG) : DENSITY
!     STZ     (NFLEVG) : GEOPOTENTIAL

! === OTHER QUANTITIES ===

! NLEXTRAP: LEVEL (CLOSE TO THE SURFACE) 
!          FROM WHERE THE TEMPERATURE IS EXTRAPOLATED 
! HEXTRAP : height from where the temperature is extrapolated.

REAL(KIND=JPRB) :: RDTDZ1
REAL(KIND=JPRB) :: RDTDZ2
REAL(KIND=JPRB) :: RDTDZ3
REAL(KIND=JPRB) :: RDTDZ4
REAL(KIND=JPRB) :: RDTDZ5
REAL(KIND=JPRB) :: RDTDZ6
REAL(KIND=JPRB) :: RDTDZ7
REAL(KIND=JPRB) :: RDTDZ8
REAL(KIND=JPRB) :: RDTDZ9
REAL(KIND=JPRB) :: RZTROP
REAL(KIND=JPRB) :: RZSTRA
REAL(KIND=JPRB) :: RZSTR2
REAL(KIND=JPRB) :: RZSTPO
REAL(KIND=JPRB) :: RZMESO
REAL(KIND=JPRB) :: RZMES2
REAL(KIND=JPRB) :: RZMEPO
REAL(KIND=JPRB) :: RZABOV
REAL(KIND=JPRB) :: RTSUR
REAL(KIND=JPRB) :: RTTROP
REAL(KIND=JPRB) :: RTSTRA
REAL(KIND=JPRB) :: RTSTR2
REAL(KIND=JPRB) :: RTSTPO
REAL(KIND=JPRB) :: RTMESO
REAL(KIND=JPRB) :: RTMES2
REAL(KIND=JPRB) :: RTMEPO
REAL(KIND=JPRB) :: RTABOV
REAL(KIND=JPRB) :: RPTROP
REAL(KIND=JPRB) :: RPSTRA
REAL(KIND=JPRB) :: RPSTR2
REAL(KIND=JPRB) :: RPSTPO
REAL(KIND=JPRB) :: RPMESO
REAL(KIND=JPRB) :: RPMES2
REAL(KIND=JPRB) :: RPMEPO
REAL(KIND=JPRB) :: RPABOV
REAL(KIND=JPRB) :: HEXTRAP

INTEGER(KIND=JPIM) :: NLEXTRAP

REAL(KIND=JPRB) :: VZTROP
REAL(KIND=JPRB) :: VZSTRA
REAL(KIND=JPRB) :: VZSTR2
REAL(KIND=JPRB) :: VZSTPO
REAL(KIND=JPRB) :: VZMESO
REAL(KIND=JPRB) :: VZMES2
REAL(KIND=JPRB) :: VZMEPO
REAL(KIND=JPRB) :: VZABOV
REAL(KIND=JPRB) :: VTSUR
REAL(KIND=JPRB) :: VTTROP
REAL(KIND=JPRB) :: VTSTRA
REAL(KIND=JPRB) :: VTSTR2
REAL(KIND=JPRB) :: VTSTPO
REAL(KIND=JPRB) :: VTMESO
REAL(KIND=JPRB) :: VTMES2
REAL(KIND=JPRB) :: VTMEPO
REAL(KIND=JPRB) :: VTABOV
REAL(KIND=JPRB) :: VPTROP
REAL(KIND=JPRB) :: VPSTRA
REAL(KIND=JPRB) :: VPSTR2
REAL(KIND=JPRB) :: VPSTPO
REAL(KIND=JPRB) :: VPMESO
REAL(KIND=JPRB) :: VPMES2
REAL(KIND=JPRB) :: VPMEPO
REAL(KIND=JPRB) :: VPABOV
REAL(KIND=JPRB) :: VDTDZ1
REAL(KIND=JPRB) :: VDTDZ2
REAL(KIND=JPRB) :: VDTDZ3
REAL(KIND=JPRB) :: VDTDZ4
REAL(KIND=JPRB) :: VDTDZ5
REAL(KIND=JPRB) :: VDTDZ6
REAL(KIND=JPRB) :: VDTDZ7
REAL(KIND=JPRB) :: VDTDZ8
REAL(KIND=JPRB) :: VDTDZ9

TYPE TSTA
REAL(KIND=JPRB), ALLOCATABLE :: STPREH(:)
REAL(KIND=JPRB), ALLOCATABLE :: STPRE(:)
REAL(KIND=JPRB), ALLOCATABLE :: STPHI(:)
REAL(KIND=JPRB), ALLOCATABLE :: STTEM(:)
REAL(KIND=JPRB), ALLOCATABLE :: STDEN(:)
REAL(KIND=JPRB), ALLOCATABLE :: STZ(:)
END TYPE TSTA

END MODULE YOMSTA
