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

MODULE YOESRTA26

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA26* - SRTM COEFFICIENTS FOR INTERVAL 26
!     BAND 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG26 = 16

REAL(KIND=JPRB) :: SFLUXREF(JPG), RAYL(JPG)

REAL(KIND=JPRB) :: SFLUXREFC(NG26), RAYLC(NG26)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
END MODULE YOESRTA26

