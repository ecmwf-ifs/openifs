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

MODULE FITSPECTRUM_MOD

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

PRIVATE
PUBLIC FITSPECTRUM, LFORCEZ, LFORCELNSP

LOGICAL :: LFORCEZ
LOGICAL :: LFORCELNSP

CONTAINS

SUBROUTINE FITSPECTRUM(KSMAX,KSPEC2G,PSPEC,KSTRUNK,PBETA0,PBETA1)



!**** *FITSPECTRUM*  - linear fit of spectrum to (n*(n+1))**P

!     Purpose.  GRIB allows for "complex packing" of spectral
!     --------  coeficients. This means that a sub-spectrum
!               is not packed and for the rest a multiplying
!               factor of (n*(n+1))**p, where n is the total
!               wave-number, is applied before packing. This routine
!               produces an optimal P (-PBETA1) based on the field
!               itself.

!**   Interface.
!     ----------
!        *CALL* *FITSPECTRUM(...)

!        Explicit arguments :   KSMAX   - spectral field truncation
!        --------------------   KSPEC2G - size of spectral field array
!                               PSPEC   - spectral field
!                               KSTRUNK - truncation of "non-packed" 
!                                         part of spectrum
!                               PBETA0  - linear fit par.  (output)
!                               PBETA1  - linear fit par.  (output)

!        Implicit arguments :  
!        --------------------

!     Method.
!     -------
!        If F(n,m) = (n*(n+1))**(-P) * G(n,m) where F(n,m) is the
!     original spectral field and n the total wavenumber the aim
!     is to minimize G in a 1 norm with respect to P.This can only
!     partially be achieved. What we do is to to compute H(n) where
!     H(n)=max(F(n,m)) with respect to m. We then perform a least 
!     square fit for the equation log(H(n))=beta0+beta1*log(n*(n+1)).
!     To ensure a better fit for the lower end of the spectrum we
!     apply an (arbitrary) weighting funtion W(n)=1./(n-kstrunk+1)
!     before fitting.

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Seber, G.A.F. (1979). Linear Regression Analyses.
!        John Wiley and Sons

!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!       Original : 94-05-05
!       Modified 95-06-06 by L.Isaksen : Reordering of spectral arrays
!       Modified 02-04-03 Y. Tremolet  : Resolution different from NSMAX
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------


INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2G 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPEC(KSPEC2G) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTRUNK 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBETA0 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBETA1 

REAL(KIND=JPRD) :: ZW(0:KSMAX),ZNORM(0:KSMAX)

INTEGER(KIND=JPIM) :: INM, ISMIN, JIR, JM, JN
INTEGER(KIND=JPIM) :: IASM0G(0:KSMAX)

REAL(KIND=JPRD) :: Z, ZEPS, ZSUM1, ZSUM2, ZWSUM, ZX, ZXMW, ZY, ZYMW, ZZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FITSPECTRUM',0,ZHOOK_HANDLE)
INM=1
DO JM=0,KSMAX
  IASM0G(JM)=INM
  INM=INM+(KSMAX+1-JM)*2
ENDDO

!*       1.    BEST FIT.
!              ---------

!*       1.1   INITIAL COMPUTATIONS

ISMIN = KSTRUNK+1
ZEPS  = EPSILON(Z)

!*       1.2   SET UP WEIGHTS

ZZ=REAL(KSMAX-ISMIN+1,JPRD)
DO JN=ISMIN,KSMAX
  ZW(JN) = ZZ/REAL(JN-ISMIN+1,JPRD)
ENDDO

!*       1.3   COMPUTE NORMS

DO JN=ISMIN,KSMAX
  ZNORM(JN) = 0.0_JPRD
ENDDO
DO JM=0,KSMAX
  DO JIR=0,MIN(1,JM)
    DO JN=MAX(ISMIN,JM),KSMAX
      INM=IASM0G(JM)+(JN-JM)*2+JIR
      ZNORM(JN) = MAX(ZNORM(JN),ABS(REAL(PSPEC(INM),JPRD)))
    ENDDO
  ENDDO
ENDDO
DO JN=ISMIN,KSMAX
  ZNORM(JN)=MAX(ZNORM(JN),ZEPS)
  IF(ZNORM(JN) == ZEPS) THEN
    ZW(JN) = 100._JPRD*ZEPS
  ENDIF
ENDDO

!*       1.4   LINEAR FIT

ZXMW  = 0.0_JPRD
ZYMW  = 0.0_JPRD
ZWSUM = 0.0_JPRD
DO JN=ISMIN,KSMAX
  ZX    = LOG(REAL(JN*(JN+1),JPRD))
  ZY    = LOG(ZNORM(JN))
  ZXMW  = ZXMW+ZX*ZW(JN)
  ZYMW  = ZYMW+ZY*ZW(JN)
  ZWSUM = ZWSUM+ZW(JN)
ENDDO
ZXMW  = ZXMW/ZWSUM
ZYMW  = ZYMW/ZWSUM
ZSUM1 = 0.0_JPRD
ZSUM2 = 0.0_JPRD
DO JN=ISMIN,KSMAX
  ZX    = LOG(REAL(JN*(JN+1),JPRD))
  ZY    = LOG(ZNORM(JN))
  ZSUM1 = ZSUM1+ZW(JN)*(ZY-ZYMW)*(ZX-ZXMW)
  ZSUM2 = ZSUM2+ZW(JN)*(ZX-ZXMW)**2
ENDDO

PBETA1 = REAL(ZSUM1/ZSUM2,JPRB)
PBETA0 = REAL(ZYMW-PBETA1*ZXMW,JPRB)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FITSPECTRUM',1,ZHOOK_HANDLE)
END SUBROUTINE FITSPECTRUM

END MODULE FITSPECTRUM_MOD
