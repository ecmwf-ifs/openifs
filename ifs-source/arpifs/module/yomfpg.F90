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

MODULE YOMFPG

USE PARKIND1, ONLY : JPIM, JPRB
USE PARFPOS,  ONLY : JPOSDOM, JPOSGL, JPOSLE

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Defining the output possibly transformed geometry
!     When not precised variables are DM-global for distributed memory.
!     Abbreviation DM-ARR stands for "DM distribution for arrival geometry"
!     Abbreviation DM-DEP stands for "DM distribution for departure geometry"

!========== FOR ALL KINDS OF OUTPUT (SUB)DOMAINS =========

! NFPMAX      : (global) truncation order
! NFPLEV      : (global) number of levels
! FPVALH      : (global) "A" coefficients of vertical system
! FPVBH       : (global) "B" coefficients of vertical system
! NFPDISTRIB  : (global) : additional transposition between departure geometry and arrival geometry in FULL-POS.
!             : = 0 : no transposition
!             : = 1 : basic transposition for gridoint outputs only
!             : = 2 : advanced transposition for compatibility with spectral transforms
! NFPFFTW     : = 0 if FFT992 should be used for the corresponding subdomain
!               = 1 if FFTW   should be used for the corresponding subdomain
! NFPFLT      : = 0 if traditional Legendre Transforms should be used for the corresponding subdomain
!               = 1 if Fast Legendre Transforms should be used for the corresponding subdomain

!========== FOR GAUSSIAN GRID AS OUTPUT ONLY =============

! FPMUCEN : MU OF THE POLE OF INTEREST
! FPLOCEN : LONGITUDE OF THE POLE OF INTEREST
! NFPHTYP : 0 = regular grid
!         : 1 = number of points proportional to sqrt(1-mu**2)
!         : 2 = number of points read on namelist namfpg
! FPNLGINC: increment to get non-linear grid
! NFPRGRI : nunber of active points on a parallel
! FPSTRET : STRETCHING FACTOR
! NFPTTYP : 1 = POLE OF STRETCHING, POLE OF THE COLLOCATION GRID
!             AT THE NORTHERN POLE OF THE REAL EARTH.
!           2 = POLE OF STRETCHING, POLE OF THE COLLOCATION GRID
!             ANYWHERE ON THE REAL EARTH.

!========== FOR ALADIN SUBDOMAIN AS OUTPUT ONLY =============

! NMFPMAX: Truncation in y
! LFPMAP : .T./.F. if the domain is defined through its coordinates/wavelengths
! FPLON0 : LA0 geographic longitude of ref. for the projection !
! FPLAT0 : FI0 geographic latitude of ref. for the projection  !  IN DEGREES !!!
! LFPMRT : If .T. in Mercator case => use of Rotated/Tilted option

TYPE TNAMFPG

INTEGER(KIND=JPIM) :: NFPMAX(JPOSDOM)
INTEGER(KIND=JPIM) :: NFPFFTW(JPOSDOM)
INTEGER(KIND=JPIM) :: NFPFLT(JPOSDOM)
INTEGER(KIND=JPIM) :: NFPRGRI(JPOSGL)
INTEGER(KIND=JPIM) :: NFPHTYP
INTEGER(KIND=JPIM) :: NFPTTYP
INTEGER(KIND=JPIM) :: NMFPMAX
REAL(KIND=JPRB) :: FPMUCEN
REAL(KIND=JPRB) :: FPLOCEN
REAL(KIND=JPRB) :: FPSTRET
REAL(KIND=JPRB) :: FPLON0
REAL(KIND=JPRB) :: FPLAT0
REAL(KIND=JPRB) :: FPNLGINC
LOGICAL :: LFPMAP
LOGICAL :: LFPMRT
INTEGER(KIND=JPIM) :: NFPDISTRIB

END TYPE TNAMFPG

TYPE TNAMFPV

INTEGER(KIND=JPIM) :: NFPLEV
REAL(KIND=JPRB) :: FPVALH(0:JPOSLE)
REAL(KIND=JPRB) :: FPVBH(0:JPOSLE)

END TYPE TNAMFPV
!     ------------------------------------------------------------------
END MODULE YOMFPG
