! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE ROTSPEC (NFRE, NANG, ML, KL, F1, F3, RTHET)

!----------------------------------------------------------------------

!**** *ROTSPEC* - ROUTINE TO ROTATE THE SPECTRUM.

!     EVA BAUER      MPI  HAMBURG    MAY 1990.
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4

!*    PURPOSE.
!     -------

!       TO ROTATE THE SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *ROTSPEC (NFRE, NANG, ML, KL, F1, F3, RTHET)*
!          *NFRE*  - FREQUENCY DIMENSION OF SPECTRA.
!          *NANG*  - DIRECTION DIMENSION OF SPECTRA.
!          *ML*    - NUMBER OF FREQUENCIES.
!          *KL*    - NUMBER OF DIRECTIONS.
!          *F1*   - SPECTRUM TO BE ROTATED.
!          *F3*   - ROTATED SPECTRUM.
!          *RTHET* - TURNING ANGLE IN RADIANS, CLOCKWISE.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ---------

!       NONE.

!     REFERENCES.
!     -----------

!       NONE.

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG, ML, KL
      REAL(KIND=JWRB), INTENT(IN) :: RTHET
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: F1
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(OUT) :: F3


      INTEGER(KIND=JWIM) :: K, M
      INTEGER(KIND=JWIM) :: INC, KC, KC1
      REAL(KIND=JWRB) :: ZPI, FTH, ADIF, BDIF

! --------------------------------------------------------------------

      ZPI=8.0_JWRB*ATAN(1._JWRB)

      FTH = MOD(RTHET+ZPI,ZPI)
      FTH = FTH * REAL(KL) / ZPI
      INC = INT(FTH)
      ADIF = FTH - INC
      BDIF = 1.0_JWRB - ADIF

      DO K=1,KL
        KC  = K  - INC
        IF (KC .LT. 1) KC  = KC + KL
        KC1 = KC - 1
        IF (KC1.LT. 1) KC1 = KC1+ KL

        DO M=1,ML
          F3(K,M) = BDIF * F1(KC,M) + ADIF * F1(KC1,M)
        ENDDO
      ENDDO

      END SUBROUTINE ROTSPEC
