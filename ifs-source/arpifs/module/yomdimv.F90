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

MODULE YOMDIMV

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

TYPE :: TDIMV

!*    Dimensions of model working arrays

! === VERTICAL RESOLUTION =====================================================

! NFLEVG : number of levels in grid point space
! NFLEVL : number of levels in Fourier and Legendre space
! NFLEVLMX : maximum NFLEVL among all PEs
! NFLSUR : over dimensioning of NFLEVL for technical reasons, always odd
! NFLSUL : number of additional levels for semi-lagrangian
! NFLSA  = 1    -NFLSUL
! NFLEN  = NFLEVG+NFLSUL
! NIOLEVG : number of levels in the whole atmosphere (used for I/Os and definitions) ; 
!           NFLEVG can be a truncation of NIOLEVG

INTEGER(KIND=JPIM) :: NFLEVG
INTEGER(KIND=JPIM) :: NFLEVL
INTEGER(KIND=JPIM) :: NFLEVLMX
INTEGER(KIND=JPIM) :: NFLSUR
INTEGER(KIND=JPIM) :: NFLSUL
INTEGER(KIND=JPIM) :: NFLSA
INTEGER(KIND=JPIM) :: NFLEN
INTEGER(KIND=JPIM) :: NIOLEVG

END TYPE TDIMV

!     ------------------------------------------------------------------

END MODULE YOMDIMV
