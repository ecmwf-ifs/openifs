! (C) Copyright 1992- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE GAWL_MOD
CONTAINS
SUBROUTINE GAWL(PFN,PL,PW,PEPS,KN,KITER,PMOD)

!**** *GAWL * - Routine to perform the Newton loop

!     Purpose.
!     --------
!           Find 0 of Legendre polynomial with Newton loop
!**   Interface.
!     ----------
!        *CALL* *GAWL(PFN,PL,PW,PEPS,KN,KITER,PMOD)

!        Explicit arguments :
!        --------------------
! PFN    Fourier coefficients of series expansion
!        for the ordinary Legendre polynomials     (in)
! PL     Gaussian latitude                         (inout)
! PW     Gaussian weight                           (out)
! PEPS   0 of the machine                          (in)
! KN     Truncation                                (in)
! KITER  Number of iterations                      (out)
! PMOD   Last modification                         (inout)

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        Newton Loop.

!     Externals.
!     ----------
!        CPLEDN

!     Reference.
!     ----------

!     ARPEGE Documentation vol.2, ch3.

!     Author.
!     -------
!        Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-12-18
!        K. Yessad (Sep 2008): cleaning, improve comments.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM

USE CPLEDN_MOD      ,ONLY : CPLEDN

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KN
REAL(KIND=JPRD),INTENT(IN)     :: PFN(0:KN/2)
REAL(KIND=JPRD),INTENT(INOUT)  :: PL
REAL(KIND=JPRD),INTENT(OUT)    :: PW
REAL(KIND=JPRD),INTENT(IN)     :: PEPS
INTEGER(KIND=JPIM),INTENT(OUT) :: KITER
REAL(KIND=JPRD),INTENT(INOUT)  :: PMOD

!     ------------------------------------------------------------------


INTEGER(KIND=JPIM) :: IFLAG, ITEMAX, JTER, IODD
REAL(KIND=JPRD) :: ZW, ZX, ZXN

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

ITEMAX = 20
ZX = PL
IFLAG = 0
IODD=MOD(KN,2)

!     ------------------------------------------------------------------

!*       2. Newton iteration.
!           -----------------

DO JTER=1,ITEMAX+1
  KITER = JTER
  CALL CPLEDN(KN,IODD,PFN,ZX,IFLAG,ZW,ZXN,PMOD)
  ZX = ZXN

  IF(IFLAG == 1) EXIT
  IF(ABS(PMOD) <= PEPS*1000._JPRD) IFLAG = 1
ENDDO

PL = ZXN
PW = ZW

!     ------------------------------------------------------------------

END SUBROUTINE GAWL
END MODULE GAWL_MOD


