! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SPR (NANG, THETAQ, THETA, ST)

! ----------------------------------------------------------------------

!**** *SPR* - ROUTINE TO COMPUTE SPREADING FACTOR.

!     SUSANNE HASSELMANN  JULY 1986.

!*    PURPOSE.
!     --------

!       COMPUTATION OF COS**2 SPREADING FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SPR (NANG, THETAQ, THETA, ST)*
!          *NANG*    INTEGER  DIMENSION FOR ANGULAR INTERVALS
!          *THETAQ*  REAL     MEAN WAVE DIRECTION
!          *THETA*   REAL     ANGULAR DIRECTIONS
!          *ST*      REAL     SPREADING FUNCTION

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.


!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : PI

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NANG
      REAL(KIND=JWRB), INTENT(IN) :: THETAQ
      REAL(KIND=JWRB), DIMENSION(NANG), INTENT(IN) :: THETA
      REAL(KIND=JWRB), DIMENSION(NANG), INTENT(OUT) :: ST

      INTEGER(KIND=JWIM) :: K

      REAL(KIND=JWRB) :: ZDP, THE

! ----------------------------------------------------------------------

      ZDP=2.0_JWRB/PI

!     SPREAD FCT. WITH COS**2.

      DO K=1,NANG
        THE = COS(THETA(K)-THETAQ)
        IF (THE > 0.0_JWRB) THEN
          ST(K) = ZDP*THE**2
          IF(ST(K).LT.0.1E-08_JWRB) ST(K)=0.0_JWRB
        ELSE
          ST(K) = 0.0_JWRB
        ENDIF
      ENDDO

      END SUBROUTINE SPR
