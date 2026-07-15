! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
! 
! THE VARIABLES X AND Y ARE THE REAL AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      USE PARKIND_WAVE, ONLY : JWIM, JWRU

      IMPLICIT NONE
      REAL(KIND=JWRU) :: X, Y, X2, Y2, RE0, IM0, RE1, IM1,              &
     &                   R1, R2, T1, T2, P1, P2, RTERM, ITERM, L

      REAL(KIND=JWRU), PARAMETER, DIMENSION(8) :: EXSQ =                &
     &   (/ 0.5641003087264E0_JWRU,  0.4120286874989E0_JWRU,            &
     &      0.1584889157959E0_JWRU,  0.3078003387255E-1_JWRU,           &
     &      0.2778068842913E-2_JWRU, 0.1000044412325E-3_JWRU,           &
     &      0.1059115547711E-5_JWRU, 0.1522475804254E-8_JWRU /)

      REAL(KIND=JWRU), PARAMETER, DIMENSION(8) :: TSQ =                 &
     &   (/ 0.0E0_JWRU,              3.19303633920635E-1_JWRU,          &
     &      1.29075862295915E0_JWRU, 2.95837445869665E0_JWRU,           &
     &      5.40903159724444E0_JWRU, 8.80407957805676E0_JWRU,           &
     &      1.34685357432515E1_JWRU, 2.02499163658709E1_JWRU /)

      INTEGER(KIND=JWIM) :: N, M, K, LL
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96E2_JWRU) GO TO 50
      IF (R2.GE.1.849E1_JWRU) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/2.0E0_JWRU
      Y2 = Y/2.0E0_JWRU
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(LOG(P1+P2)/2.0E0_JWRU+0.5772156649015329E0_JWRU)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -ATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0E0_JWRU
      ITERM = 0.0E0_JWRU
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5E0_JWRU
      RE1 = T1
      IM1 = T2
      P2 = SQRT(R2)
      L = 2.106E0_JWRU*P2 + 4.4E0_JWRU
      IF (P2.LT.8.0E-1_JWRU) L = 2.129E0_JWRU*P2 + 4.0E0_JWRU
      LL=NINT(L)
      DO N=1,LL
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5E0_JWRU/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0E0_JWRU
        T1 = T1 + 0.5E0_JWRU/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
      END DO
      R1 = X/R2 - 0.5E0_JWRU*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5E0_JWRU*(X*IM1+Y*RE1)
      P1 = EXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = 2.0E0_JWRU*X
      Y2 = 2.0E0_JWRU*Y
      R1 = Y2*Y2
      P1 = SQRT(X2*X2+R1)
      P2 = SQRT(P1+X2)
      T1 = EXSQ(1)/(2.0E0_JWRU*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0E0_JWRU
      IM1 = 0.0E0_JWRU
      DO N=2,8
        T2 = X2 + TSQ(N)
        P1 = SQRT(T2*T2+R1)
        P2 = SQRT(P1+T2)
        T1 = EXSQ(N)/P1
        RE0 = RE0 + T1*P2
        IM0 = IM0 + T1/P2
        T1 = EXSQ(N)*TSQ(N)
        RE1 = RE1 + T1*P2
        IM1 = IM1 + T1/P2
      END DO
      T2 = -Y2*IM0
      RE1 = RE1/R2
      R2 = Y2*IM1/R2
      RTERM = 1.41421356237309E0_JWRU*COS(Y)
      ITERM = -1.41421356237309E0_JWRU*SIN(Y)
! THE CONSTANT IN THE PREVIOUS STATEMENTS IS,OF COURSE,
! SQRT(2.0).
      IM0 = RE0*ITERM + T2*RTERM
      RE0 = RE0*RTERM - T2*ITERM
      T1 = RE1*RTERM - R2*ITERM
      T2 = RE1*ITERM + R2*RTERM
      RE1 = T1*X + T2*Y
      IM1 = -T1*Y + T2*X
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE
! ASYMPTOTIC EXPANSIONS
   50 RTERM = 1.0E0_JWRU
      ITERM = 0.0E0_JWRU
      RE0 = 1.0E0_JWRU
      IM0 = 0.0E0_JWRU
      RE1 = 1.0E0_JWRU
      IM1 = 0.0E0_JWRU
      P1 = 8.0E0_JWRU*R2
      P2 = SQRT(R2)
      L = 3.91E0_JWRU+8.12E1_JWRU/P2
      LL=NINT(L)
      R1 = 1.0E0_JWRU
      R2 = 1.0E0_JWRU
      M = -8
      K = 3
      DO N=1,LL
        M = M + 8
        K = K - M
        R1 = REAL(K-4,JWRU)*R1
        R2 = REAL(K,JWRU)*R2
        T1 = REAL(N,JWRU)*P1
        T2 = RTERM
        RTERM = (T2*X+ITERM*Y)/T1
        ITERM = (-T2*Y+ITERM*X)/T1
        RE0 = RE0 + R1*RTERM
        IM0 = IM0 + R1*ITERM
        RE1 = RE1 + R2*RTERM
        IM1 = IM1 + R2*ITERM
      END DO
      T1 = SQRT(P2+X)
      T2 = -Y/T1
      P1 = 8.86226925452758E-1_JWRU/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*COS(Y)
      ITERM = -P1*SIN(Y)
      R1 = RE0*RTERM - IM0*ITERM
      R2 = RE0*ITERM + IM0*RTERM
      RE0 = T1*R1 - T2*R2
      IM0 = T1*R2 + T2*R1
      R1 = RE1*RTERM - IM1*ITERM
      R2 = RE1*ITERM + IM1*RTERM
      RE1 = T1*R1 - T2*R2
      IM1 = T1*R2 + T2*R1

      END SUBROUTINE KZEONE

