! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

REAL(KIND=JWRB) FUNCTION W_MODE_ST (RN3, RN2, RN1)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     DETERMINE MODE OF THE SPACE-TIME PROBABILITY DISTRIBUTION
!                                                                              !                                                                              
!     PURPOSE.
!     --------
!                                                                              !
!           DETERMINATION OF THE MODE OF THE SPACE-TIME PROBABILITY DISTR.     !
!           USING THE NEWTON METHOD

!           SOLVE FOR X: (R3*X**2+R2*X+R1)*EXP(-0.5*X**2)=1

!           Barbariol et al., 2019:
!           Maximum wave height from global model reanalysis
!                                                                              !
! ---------------------------------------------------------------------------- !

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

IMPLICIT NONE

!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
!                                                                              !
REAL(KIND=JWRB) :: RN3   !! AVERAGE NUMBER OF 3D WAVES
REAL(KIND=JWRB) :: RN2   !! AVERAGE NUMBER OF 2D WAVES
REAL(KIND=JWRB) :: RN1   !! AVERAGE NUMBER OF 1D WAVES


! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER(KIND=JWIM), PARAMETER :: MITER = 20        ! MAX NUMBER OF ITERATIONS
INTEGER(KIND=JWIM) :: ITER

REAL(KIND=JWRB) :: Z0, RES, FPRIME
REAL(KIND=JWRB), PARAMETER :: TOL = 1.0E-6_JWRB    ! TOLERANCE

!     INLINE FUNCTIONS
!     ----------------
REAL(KIND=JWRB) :: F, DF, X, R1, R2, R3
! Probabilty of maximum crest height is equal to 1
F(X,R1,R2,R3) = (X*(R3*X+R2)+R1)*EXP(-0.5_JWRB*X**2)-1._JWRB
! Derivative of F:
DF(X,R1,R2,R3) = (-X**2*(R3*X+R2)+(2_JWRB*R3-R1)*X+R1+R2)*EXP(-0.5_JWRB*X**2)


! ---------------------------------------------------------------------------- !


! FIRST GUESS: MODE OF THE TIME PROBABIILITY DISTRIBUTION
!Z0 = SQRT(2._JWRB*LOG(RN1))
Z0 = SQRT(2._JWRB*LOG(RN3)+2._JWRB*LOG(2._JWRB*LOG(RN3)+2._JWRB*LOG(2._JWRB*LOG(RN3))))
ITER = 0
RES = ABS(F(Z0, RN1, RN2, RN3))

DO WHILE (ITER < MITER .AND. TOL < RES)
   ITER = ITER + 1
   FPRIME = DF(Z0, RN1, RN2, RN3)
   IF(FPRIME /= 0._JWRB) Z0 = Z0 - F(Z0, RN1, RN2, RN3)/FPRIME
   RES = ABS(F(Z0, RN1, RN2, RN3))
ENDDO
W_MODE_ST = Z0


END FUNCTION W_MODE_ST
