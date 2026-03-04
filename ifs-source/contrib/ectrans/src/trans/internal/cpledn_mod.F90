! (C) Copyright 1987- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE CPLEDN_MOD
CONTAINS
SUBROUTINE CPLEDN(KN,KODD,PFN,PX,KFLAG,PW,PXN,PXMOD)

!**** *CPLEDN* - Routine to perform a single Newton iteration step to find
!                the zero of the ordinary Legendre polynomial of degree N

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *CPLEDN(KN,KDBLE,PX,KFLAG,PW,PXN,PXMOD)*

!        Explicit arguments :
!        --------------------
!          KN       :  Degree of the Legendre polynomial              (in)
!          KODD     :  odd or even number of latitudes                (in)
!          PFN      :  Fourier coefficients of series expansion       (in)
!                      for the ordinary Legendre polynomials
!          PX       :  abcissa where the computations are performed   (in)
!          KFLAG    :  When KFLAG.EQ.1 computes the weights           (in)
!          PW       :  Weight of the quadrature at PXN                (out)
!          PXN      :  new abscissa (Newton iteration)                (out)
!          PXMOD    :  PXN-PX                                         (out)

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        Michel Rochas, 90-08-30 (Lobatto+cleaning)
!        K. Yessad (Sep 2008): cleaning, improve comments.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KN
INTEGER(KIND=JPIM),INTENT(IN) :: KODD
REAL(KIND=JPRD),INTENT(IN)    :: PFN(0:KN/2)
REAL(KIND=JPRD),INTENT(IN)    :: PX
INTEGER(KIND=JPIM),INTENT(IN) :: KFLAG
REAL(KIND=JPRD),INTENT(OUT)   :: PW
REAL(KIND=JPRD),INTENT(INOUT) :: PXN
REAL(KIND=JPRD),INTENT(OUT)   :: PXMOD

!     ------------------------------------------------------------------

REAL(KIND=JPRD) :: ZDLX,ZDLK,ZDLLDN,ZDLXN,ZDLMOD

INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(PX)

INTEGER(KIND=JPIM) :: JN, IK

!      -----------------------------------------------------------------

!*       1. NEWTON ITERATION STEP.
!           ----------------------

ZDLX = PX

ZDLK = 0.0_JPRD
IF( KODD==0 ) ZDLK=0.5_JPRD*PFN(0)
ZDLXN = 0.0_JPRD
ZDLLDN = 0.0_JPRD
IK=1

IF(KFLAG == 0)THEN
  DO JN=2-KODD,KN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(IK)*COS(REAL(JN,JPKD)*ZDLX)
    ! normalised derivative == d/d\theta(\overbar{P_n}^0)
    ZDLLDN = ZDLLDN - PFN(IK)*REAL(JN,JPKD)*SIN(REAL(JN,JPKD)*ZDLX)
    IK=IK+1
  ENDDO
  ! Newton method
  ZDLMOD = -ZDLK/ZDLLDN
  ZDLXN = ZDLX+ZDLMOD
  PXN = ZDLXN
  PXMOD = ZDLMOD
ENDIF

!     ------------------------------------------------------------------

!*         2.    Computes weight.
!                ----------------

IF(KFLAG == 1)THEN
  DO JN=2-KODD,KN,2
    ! normalised derivative
    ZDLLDN = ZDLLDN - PFN(IK)*REAL(JN,JPKD)*SIN(REAL(JN,JPKD)*ZDLX)
    IK=IK+1
  ENDDO
  PW = REAL(2*KN+1,JPKD)/ZDLLDN**2
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE CPLEDN
END MODULE CPLEDN_MOD
