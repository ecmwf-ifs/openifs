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

SUBROUTINE SUVFE_CPSPLINES(KORDER,KBASIS,PTM,PKNOT,PSPLINE,PETAMAX)

!**** *SUVFE_CPSPLINES*  - compute coefficients
!                          of splines of given order above given
!                          sequence of knots.
!                          Used in setup of VFE scheme.

!**   Interface.
!     ----------

!     *CALL* SUVFE_CPSPLINES

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        KORDER  - order of splines 
!        KBASIS  - number of splines to be computed
!        PTM     - transformation matrix (from sampled values to local coordinate)
!        PKNOT   - sequence of nondecreasing data with dimention KORDER+KBASIS
!                  * splines consists of piecewise polynomial of order KORDER 
!                    ( linear spline has order 2 and cubis one has order 4)
!                    and these polynomials connect at PKNOT points
!                  * PKNOT are expressed in global coordinate eta

!      * OUTPUT:
!        PSPLINE - coefficients of computed splines on intervals <knot(i),knot(i+1)>;
!                  the corresponding value of i-th spline and its j-th piecewise polynomial 
!                  is S(i,j) = SUM(K,1,KORDER) PSPLINE(i,j,K)*t^(K-1)
!        PETAMAX - position of maxima of splines in eta coordinate

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Jozef Vivoda SHMU/LACE

!     Modifications.
!     --------------
!      Original : 2009-10
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      P.Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE SUVFE_HLP ,ONLY : SAMPLE_VALUE, FT2X

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KORDER
INTEGER(KIND=JPIM),INTENT(IN) :: KBASIS
REAL(KIND=JPRB),INTENT(IN)    :: PTM    (KORDER,KORDER)
REAL(KIND=JPRB),INTENT(IN)    :: PKNOT  (KBASIS + KORDER)
REAL(KIND=JPRB),INTENT(OUT)   :: PSPLINE(KBASIS,KORDER,KORDER)
REAL(KIND=JPRB),INTENT(OUT)   :: PETAMAX(KBASIS)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IKNOT, JI, ISEG, JK
REAL(KIND=JPRB)    :: ZT(KORDER)
REAL(KIND=JPRB)    :: ZSAMPLE(KORDER,KBASIS)
REAL(KIND=JPRB)    :: ZD, ZETA, ZLIM, ZMINDETA
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suvfe_basis.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_CPSPLINES',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

! amount of knots
IKNOT = KBASIS + KORDER

! limit for maxima search
ZLIM = 1.0E-4

! set all polynomials to 0.0
PSPLINE = 0.0_JPRB

! ZT = uniformly sampled local coordinate on interval <0,1>
DO JI=1,KORDER
  ZT(JI) = SAMPLE_VALUE(KORDER,JI)
ENDDO

ZMINDETA = 0.0_JPRB

! compute polynomial coefficients above each interval (i,i+1),i=1,IKNOT-1
! the polynomial coefficient are set to 0 for intervals
! with zero length (multiple nodes)
DO JK = 1, IKNOT-1  ! JK - k-th interval of domain, JK=1,IKNOT-1

  ! distance of knots at k-th interval
  ZD = PKNOT(JK+1) - PKNOT(JK)

  ! skip intervals with zero size
  IF(ABS(ZD) > ZMINDETA)THEN

    ! sample interval (JK,JK+1) with KORDER values
    ! determined by ZT array
    DO JI = 1,KORDER

      ! transform local interval coordinate ZT (t=<0,1>)
      ! into global variable eta
      ! ZETA = ZT(JI)*ZD + PKNOT(JK)
      ZETA = FT2X(PKNOT(JK), PKNOT(JK+1), ZT(JI))

      ! evalute all splines on domain in global coordinate ZETA
      CALL SUVFE_BASIS(KORDER,ZETA,KBASIS,PKNOT,ZSAMPLE(JI,:))

    ENDDO

    DO JI = MAX(1,JK-KORDER+1),MIN(KBASIS,JK)  ! JI - i-th basis spline function
      ISEG = JK + 1 - JI   ! ISEG - i-th segment of piecewise spline basis function
      PSPLINE(JI,ISEG,:) = MATMUL(PTM,ZSAMPLE(:,JI))
    ENDDO

  ENDIF

ENDDO

! greville definition of eta levels from KNOT positions
DO JI = 1, KBASIS
  PETAMAX(JI) = 0.0_JPRB
  DO JK = 1, KORDER - 1
    PETAMAX(JI) = PETAMAX(JI) + PKNOT(JI+JK)
  ENDDO
  PETAMAX(JI) = PETAMAX(JI) / REAL(KORDER - 1)
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_CPSPLINES',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_CPSPLINES
