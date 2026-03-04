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

SUBROUTINE SUVFE_CPSPLINES(YDCVER,KORDER,KBASIS,PTM,PKNOT,PSPLINE,PETAMAX)

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
!   Original : 2009-10
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMCVER,ONLY : TCVER
USE YOMLUN    ,ONLY : NULOUT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE SUVFE_HLP ,ONLY : SAMPLE_VALUE, DPOL, EVPOL, FT2X, IPOL, &
                    & RTMIN, RTMAX, SETTM, FX2T, GLOBAL2LOCAL, MULPOL

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCVER)    ,INTENT(IN)    :: YDCVER
INTEGER(KIND=JPIM),INTENT(IN) :: KORDER
INTEGER(KIND=JPIM),INTENT(IN) :: KBASIS
REAL(KIND=JPRB),INTENT(IN)    :: PTM    (KORDER,KORDER)
REAL(KIND=JPRB),INTENT(IN)    :: PKNOT  (KBASIS + KORDER)
REAL(KIND=JPRB),INTENT(OUT)   :: PSPLINE(KBASIS,KORDER,KORDER)
REAL(KIND=JPRB),INTENT(OUT)   :: PETAMAX(KBASIS)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IKNOT, JI, ISEG, JK, ITER
REAL(KIND=JPRB)    :: ZT     (KORDER)
REAL(KIND=JPRB)    :: ZSAMPLE(KORDER,KBASIS)
REAL(KIND=JPRB)    :: ZD, ZETA, ZLIM, ZX1, ZX2, ZFAC
REAL(KIND=JPRB)    :: ZPD1(KORDER-1), ZPD2(KORDER-2), ZINC
REAL(KIND=JPRB)    :: ZPI(KORDER+1), ZP(KORDER), ZPM(2*KORDER), ZPMI(2*KORDER+1)
REAL(KIND=JPRB)    :: ZD1, ZD2, ZVAL, ZZ, ZMAX(KBASIS), ZTM(KORDER,KORDER)
REAL(KIND=JPRB)    :: ZETA1, ZETA2, ZT1, ZT2, ZCONST, ZCIF
!REAL(KIND=JPRB)    :: EVPOL
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

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

ZMAX = 0.0_JPRB

! compute polynomial coefficients above each interval (i,i+1),i=1,IKNOT-1
! the polynomial coefficient are set to 0 for intervals
! with zero length (multiple nodes)
DO JK = 1, IKNOT-1  ! JK - k-th interval of domain, JK=1,IKNOT-1

  ! distance of knots at k-th interval
  ZD = PKNOT(JK+1) - PKNOT(JK)

  ! skip intervals with zero size
  IF(ABS(ZD) > YDCVER%RMINDETA)THEN

    ! sample interval (JK,JK+1) with KORDER values
    ! determined by ZT array
    DO JI = 1,KORDER

      ! transform local interval coordinate ZT (t=<0,1>)
      ! into global variable eta
      ZETA = FT2X(PKNOT(JK), PKNOT(JK+1), ZT(JI))

      ! evalute all splines on domain in global coordinate ZETA
      CALL SUVFE_BASIS(KORDER,ZETA,KBASIS,PKNOT,ZSAMPLE(JI,:))

    ENDDO

    DO JI = MAX(1,JK-KORDER+1),MIN(KBASIS,JK)  ! JI - i-th basis spline function
      ISEG = JK + 1 - JI      ! ISEG - i-th segment of piecewise spline basis function

      PSPLINE(JI,ISEG,:) = MATMUL(PTM,ZSAMPLE(:,JI))

      IF(YDCVER%LVFE_VERBOSE)THEN

        !-------------------------
        ! DBG PART

        CALL SETTM(PKNOT(JK), PKNOT(JK+1), KORDER, ZTM)
        ZP = MATMUL(ZTM,ZSAMPLE(:,JI))

        ZETA1 = PKNOT(JK) + 0.25_JPRB * ZD
        ZETA2 = PKNOT(JK) + 0.75_JPRB * ZD
        ZT1   = FX2T(PKNOT(JK), PKNOT(JK+1), ZETA1)
        ZT2   = FX2T(PKNOT(JK), PKNOT(JK+1), ZETA2)

        ! value of spline at the beginning of knot interval
        ZCONST = 0.0_JPRB
        CALL IPOL(KORDER,ZP,ZCONST,KORDER+1,ZPI)
        CALL MULPOL(KORDER, ZP, KORDER + 1, ZPI, 2*KORDER, ZPM)
        CALL IPOL(2*KORDER,ZPM,0.0_JPRB,2*KORDER+1,ZPMI)
        ZX1 = EVPOL(2*KORDER+1, ZPMI, ZETA1)
        ZX2 = EVPOL(2*KORDER+1, ZPMI, ZETA2)
        ZFAC = 1.0_JPRB
        ! WRITE(NULOUT,'("DBG ETA INTEGRAL :: ",I3,1X,I3,3(1X,F15.12))') JI, ISEG, ZFAC * (ZX2 - ZX1)

        ! value of spline at the beginning of knot interval
        CALL IPOL(KORDER,PSPLINE(JI,ISEG,:),ZCONST,KORDER+1,ZPI)

        ! this computation must be done in T coordinate of given function

        ZX1 = EVPOL(KORDER+1, ZPI, RTMIN)
        ZX2 = EVPOL(KORDER+1, ZPI, ZT1)
        ZCIF = ZD * (ZX2 - ZX1)

        ZX1 = EVPOL(KORDER+1, ZPI, ZT1)
        ZX2 = EVPOL(KORDER+1, ZPI, ZT2)
        ZCONST = ZD * (ZX2 - ZX1)

        CALL MULPOL(KORDER,PSPLINE(JI,ISEG,:),KORDER+1,ZPI,2*KORDER,ZPM)
        CALL IPOL(2*KORDER,ZPM,0.0_JPRB,2*KORDER+1,ZPMI)
        ZX1 = EVPOL(2*KORDER+1, ZPMI, ZT1)
        ZX2 = EVPOL(2*KORDER+1, ZPMI, ZT2)
        ZFAC = ZD*ZD
        ! WRITE(NULOUT,'("DBG T01 INTEGRAL :: ",I3,1X,I3,3(1X,F15.12))') JI, ISEG, ZFAC * (ZX2 - ZX1)

        ! this interval is not the same is previous two
        CALL GLOBAL2LOCAL(KORDER,PTM,ZT1,ZT2,PSPLINE(JI,ISEG,:),ZP)
        CALL IPOL(KORDER,ZP,0.0_JPRB,KORDER+1,ZPI)

        CALL MULPOL(KORDER,ZP,KORDER+1,ZPI,2*KORDER,ZPM)
        CALL IPOL(2*KORDER,ZPM,0.0_JPRB,2*KORDER+1,ZPMI)
        ZX1 = EVPOL(2*KORDER+1, ZPMI, RTMIN)
        ZX2 = EVPOL(2*KORDER+1, ZPMI, RTMAX)
        ZFAC = (ZT2 - ZT1) * (ZT2 - ZT1) * ZD * ZD

        ! WRITE(NULOUT,'("DBG LOC INTEGRAL :: ",I3,1X,I3,3(1X,F15.12))') JI, ISEG, ZFAC * (ZX2 - ZX1) + ZCONST * ZCIF, ZCONST, ZFAC

      ENDIF

      !-------------------------

      IF(YDCVER%LVFE_MAXIMAS)THEN

        !---------------
        ! compute maxima of spline of this interval
        ! newton method:
        !     t(n+1) = t(n) - f'(t(n)) / f''(t(n))
        !
        ! initial value: t(0) = 1/2
        !---------------
        CALL DPOL(KORDER  ,PSPLINE(JI,ISEG,:),KORDER-1,ZPD1 )
        CALL DPOL(KORDER-1,ZPD1              ,KORDER-2,ZPD2)
        ZVAL = 0.5_JPRB * (RTMIN + RTMAX)
        DO ITER = 1, 20
          ZD1     = EVPOL(KORDER-1, ZPD1, ZVAL) / ZD
          ZD2     = EVPOL(KORDER-2, ZPD2, ZVAL) / (ZD**2)
          ZINC    = ZD1 / ZD2
          IF(ABS(ZINC) > ZLIM)THEN
            ZVAL = ZVAL - ZINC
          ELSE
            EXIT
          ENDIF
        ENDDO

        ZZ = EVPOL(KORDER ,PSPLINE(JI,ISEG,:),    ZVAL)
        IF(ZZ > ZMAX(JI) .AND. (ZVAL > RTMIN .AND. ZVAL < RTMAX ))THEN
          PETAMAX(JI) = FT2X(PKNOT(JK), PKNOT(JK+1), ZVAL)
          ZMAX (JI) = ZZ
        ENDIF

        ! check BC at t=0.0
        ZZ = EVPOL(KORDER ,PSPLINE(JI,ISEG,:),RTMIN)
        IF(ZZ > ZMAX(JI))THEN
          PETAMAX(JI) = PKNOT(JK)
          ZMAX (JI) = ZZ
        ENDIF

        ! check BC at t=1.0
        ZZ = EVPOL(KORDER ,PSPLINE(JI,ISEG,:),RTMAX)
        IF(ZZ > ZMAX(JI))THEN
          PETAMAX(JI) = ZD + PKNOT(JK)
          ZMAX (JI) = ZZ
        ENDIF

      ENDIF

    ENDDO

  ENDIF

ENDDO


IF (.NOT.YDCVER%LVFE_MAXIMAS) THEN

  ! greville definition of eta levels
  ! from KNOT positions

  DO JI = 1, KBASIS
    PETAMAX(JI) = 0.0_JPRB
    DO JK = 1, KORDER - 1
      PETAMAX(JI) = PETAMAX(JI) + PKNOT(JI+JK)
    ENDDO
    PETAMAX(JI) = PETAMAX(JI) / REAL(KORDER - 1)
  ENDDO

ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_CPSPLINES',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_CPSPLINES
