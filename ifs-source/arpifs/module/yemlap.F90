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

MODULE YEMLAP

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    Constants related to the Laplace space

!  NCPL2M :  number of complex Laplace coefficient for m given (local)
!  NCPL4M :  number of real    Laplace coefficient for m given (local)
!  NCPLM  :  number of complex real/imag. Laplace coefficient for m given
!            (half of NCPL2M)
!  NCPL2N :  number of complex Laplace coefficient for n given
!  NCPL4N :  number of real    Laplace coefficient for n given
!  NCPLN  :  number of complex real/imag. coefficient for n given
!            (half of NCPL2N)

!  RLEPDIN:  eigen-values of the Laplace operator
!  RLEPINN:  eigen-values of its inverse
!  RLEPDIM:  eigen-values of the Laplace operator
!  RLEPINM:  eigen-values of its inverse

!  NESM0 : address (local) in a spectral array of (n, m=0)
!  NESM0G: address (global) in a spectral array of (n, m=0)
!  NESPZERO: global adresses in the (global-dimensioned) spectral arrays
!            of those coefficients which are always zero, id est:
!            m = 0 - coefficients 3 and 4 for every quadruplet (1,2,3,4)
!            n = 0 - coefficients 2 and 4 for every quadruplet (1,2,3,4)

!  NPME  : address for the Laplace operator and its inverse
!  NPNE  : address for the Laplace operator and its inverse
!  MVALUE: wavenumber M for the given spectral coefficient

TYPE :: TLEP
  INTEGER(KIND=JPIM),ALLOCATABLE:: NCPL2M(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NCPL4M(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NCPLM(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NCPL2N(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NCPL4N(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NCPLN(:)

  REAL(KIND=JPRB),ALLOCATABLE:: RLEPDIN(:)
  REAL(KIND=JPRB),ALLOCATABLE:: RLEPINN(:)
  REAL(KIND=JPRB),ALLOCATABLE:: RLEPDIM(:)
  REAL(KIND=JPRB),ALLOCATABLE:: RLEPINM(:)

  INTEGER(KIND=JPIM),ALLOCATABLE:: NESM0(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NESPZERO(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NESM0G(:)

  INTEGER(KIND=JPIM),ALLOCATABLE:: NPME(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: NPNE(:)
  INTEGER(KIND=JPIM),ALLOCATABLE:: MVALUE(:)
END TYPE TLEP

!     ------------------------------------------------------------------
END MODULE YEMLAP
