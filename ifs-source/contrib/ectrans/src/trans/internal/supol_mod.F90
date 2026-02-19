! (C) Copyright 1987- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUPOL_MOD
CONTAINS
SUBROUTINE SUPOL(KNSMAX,PDDMU,PFN,PDDPOL)

!**** *SUPOL * - Routine to compute the Legendre polynomials

!     Purpose.
!     --------
!           For a given value of mu, computes the Legendre polynomials.

!**   Interface.
!     ----------
!        *CALL* *SUPOL(...)

!        Explicit arguments :
!        --------------------
!        KNSMAX    :  Truncation  (triangular)                            [in]
!        PDDMU     :  Abscissa at which the polynomials are computed (mu) [in]
!        PFN       :  Fourier coefficients of series expansion 
!                     for the ordinary Legendre polynomials               [in]
!        PDDPOL    :  Polynomials (the first index is m and the second n) [out]

!        Implicit arguments :   None
!        --------------------

!     Method.
!     -------
!        See documentation about spectral transforms
!         (doc (IDTS) by K. Yessad, appendix 3, or doc (NTA30) by M. Rochas)

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        K. YESSAD (MAY 1998): modification to avoid underflow.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib 11-Apr-2007 Emulation of vectorized quadruple precision
!                                 on NEC
!        K. YESSAD (NOV 2008): make consistent arp/SUPOLA and tfl/SUPOL.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!        R. El Khatib 30-Apr-2013 Open-MP parallelization
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM
USE TPM_POL

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KNSMAX
REAL(KIND=JPRD)   ,INTENT(IN)  :: PDDMU
REAL(KIND=JPRD)   ,INTENT(IN)  :: PFN(0:KNSMAX,0:KNSMAX)

REAL(KIND=JPRD)   ,INTENT(OUT) :: PDDPOL(0:KNSMAX,0:KNSMAX)

REAL(KIND=JPRD) :: ZDLX,ZDLX1,ZDLSITA,ZDL1SITA,ZDLS,ZDLK,ZDLLDN

INTEGER(KIND=JPIM) :: JM, JN, JK
REAL(KIND=JPRD) :: Z
REAL(KIND=JPRD) :: DCL, DDL

!     ------------------------------------------------------------------

!*       1. First two columns.
!           ------------------

ZDLX=PDDMU
ZDLX1=ACOS(ZDLX)
ZDLSITA=SQRT(1.0_JPRD-ZDLX*ZDLX)

PDDPOL(0,0)=1._JPRD
ZDLLDN = 0.0_JPRD

! IF WE ARE LESS THAN 1Meter FROM THE POLE,
IF(ABS(REAL(ZDLSITA,KIND(Z))) <= SQRT(EPSILON(Z)))THEN
  ZDLX=1._JPRD
  ZDLSITA=0._JPRD
  ZDL1SITA=0._JPRD
ELSE
  ZDL1SITA=1.0_JPRD/ZDLSITA
ENDIF

!*          ordinary Legendre polynomials from series expansion
!           ---------------------------------------------------

! even N
!$OMP PARALLEL DO PRIVATE(JN,ZDLK,ZDLLDN,JK)
DO JN=2,KNSMAX,2
  ZDLK = 0.5_JPRD*PFN(JN,0)
  ZDLLDN = 0.0_JPRD
  ! represented by only even k
  DO JK=2,JN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(JN,JK)*COS(DDI(JK)*ZDLX1)
    ! normalised associated Legendre polynomial == \overbar{P_n}^1
    ZDLLDN = ZDLLDN + DDA(JN)*PFN(JN,JK)*DDI(JK)*SIN(DDI(JK)*ZDLX1)
  ENDDO
  PDDPOL(0,JN) = ZDLK
  PDDPOL(1,JN) = ZDLLDN
ENDDO
!$OMP END PARALLEL DO
! odd N
!$OMP PARALLEL DO PRIVATE(JN,ZDLK,ZDLLDN,JK)
DO JN=1,KNSMAX,2
  ZDLK = 0.0_JPRD
  ZDLLDN = 0.0_JPRD
  ! represented by only odd k
  DO JK=1,JN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(JN,JK)*COS(DDI(JK)*ZDLX1)
    ! normalised associated Legendre polynomial == \overbar{P_n}^1
    ZDLLDN = ZDLLDN + DDA(JN)*PFN(JN,JK)*DDI(JK)*SIN(DDI(JK)*ZDLX1)
  ENDDO
  PDDPOL(0,JN) = ZDLK
  PDDPOL(1,JN) = ZDLLDN
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

!*       2. Diagonal (the terms 0,0 and 1,1 have already been computed)
!           Belousov, equation (23)
!           -----------------------------------------------------------

ZDLS=ZDL1SITA*TINY(ZDLS)

#ifdef VPP
!OCL SCALAR
#endif
DO JN=2,KNSMAX
  PDDPOL(JN,JN)=PDDPOL(JN-1,JN-1)*ZDLSITA*DDH(JN)
  IF ( ABS(PDDPOL(JN,JN)) < ZDLS ) PDDPOL(JN,JN)=0.0_JPRD
ENDDO

!     ------------------------------------------------------------------

!*       3. General recurrence (Belousov, equation 17)
!           -----------------------------------------

DO JN=3,KNSMAX
!DIR$ IVDEP
!OCL NOVREC
  DO JM=2,JN-1
    PDDPOL(JM,JN)=DDC(JM,JN)*PDDPOL(JM-2,JN-2)&
     &-DDD(JM,JN)*PDDPOL(JM-2,JN-1)*ZDLX &
     &+DDE(JM,JN)*PDDPOL(JM  ,JN-1)*ZDLX
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SUPOL
END MODULE SUPOL_MOD
