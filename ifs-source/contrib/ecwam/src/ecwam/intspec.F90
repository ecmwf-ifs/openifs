! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INTSPEC (NFRE, NANG, ML, KL, FR, DEL12, DEL1L,         &
     &                    F1, FMEAN1, EMEAN1, THETM1,                   &
     &                    F2, FMEAN2, EMEAN2, THETM2,                   &
     &                    FL, FMEAN,  EMEAN,  THETM )

! ----------------------------------------------------------------------

!**** *INTSPEC* -  INTERPOLATION OF SPECTRA.

!     SUSANNE HASSELMANN  MPI        JUNE 1990.
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4

!*    PURPOSE.
!     --------

!       INTERPOLATION OF SPECTRA.

!**   INTERFACE.
!     ----------

!       *CALL* *INTSPEC (NFRE, NANG, ML, KL, FR, DEL12, DEL1L,
!                        F1, FMEAN1, EMEAN1, THETM1,
!                        F2, FMEAN2, EMEAN2, THETM2,
!                        FL, FMEAN,  EMEAN,  THETM )*
!         *NFRE*   - FREQUENCY DIMENSION OF SPECTRA.
!         *NANG*   - DIRECTION DIMENSION OF SPECTRA.
!         *ML*     - NUMBER OF FREQUENCIES.
!         *KL*     - NUMBER OF DIRECTIONS.
!         *FR*     - FREQUENCY ARRAY.
!         *DEL12*  - DISTANCE SPECTRUM 2 - SPECTRUM 1.
!         *DEL1L*  - DISTANCE SPECTRUM L - SPECTRUM 1.
!         *F1*     - INPUT SPECTRUM 1.
!         *FMEAN1* - INPUT MEAN FREQUENCY OF F1.
!         *EMEAN1* - INPUT MEAN ENERGY OF F1.
!         *THETM1* - INPUT MEAN DIRECTION OF F1.
!         *F2*     - INPUT SPECTRUM 2.
!         *FMEAN2* - INPUT MEAN FREQUENCY OF F2.
!         *EMEAN2* - INPUT MEAN ENERGY OF F2.
!         *THETM2* - INPUT MEAN DIRECTION OF F2.
!         *FL*     - INTEPOLATED SPECTRUM.
!         *FMEAN*  - INTEPOLATED MEAN FREQUENCY.
!         *EMEAN*  - INTEPOLATED MEAN ENERGY.
!         *THETM*  - INTEPOLATED MEAN DIRECTION.

!     METHOD.
!     -------

!       ROTATE SPECTRA ACCORDING TO MEAN OF MEAN ANGLES, TRANSFORM
!       FREQUENCIES ACCORDING TO MEAN OF MEAN FREQUENCIES ,ADJUST ENERGY
!       ACCORDCING TO MEAN OF TOTAL ENERGY AND INTERPOLATE RESULTING
!       SPECTRA.

!     EXTERNALS.
!     ----------

!       *ROTSPEC*   - ROTATES SPECTRUM.
!       *STRSPEC*   - TRANSFORM FREQUENCIES.

!     REFERENCES.
!     -----------

!       K.HASSELMANN, 1990,
!          INTERPOLATION OF WAVE SPECTRA. WAM NOTE 6/6/90.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : ZPI
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "rotspec.intfb.h"
#include "strspec.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG, ML, KL
      REAL(KIND=JWRB), INTENT(IN) :: DEL12,  DEL1L
      REAL(KIND=JWRB), INTENT(IN) :: FMEAN1, EMEAN1, THETM1
      REAL(KIND=JWRB), INTENT(IN) :: FMEAN2, EMEAN2, THETM2
      REAL(KIND=JWRB), INTENT(OUT) :: FMEAN, EMEAN, THETM
      REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(IN) :: FR
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: F1
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: F2
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(OUT) :: FL

      INTEGER(KIND=JWIM) :: M, K

      REAL(KIND=JWRB) :: GW1, GW2 
      REAL(KIND=JWRB) :: CM1, CM2, SM1, SM2, CM,  SM 
      REAL(KIND=JWRB) :: RTHET1,RTHET2, GAMMA, EMEANH
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: F3, F4

! ----------------------------------------------------------------------

!*    1. INTERPOLATION WEIGHTS.
!        ----------------------

      GW1 = (DEL12-DEL1L)/DEL12
      GW2 = DEL1L/DEL12

! ----------------------------------------------------------------------

!*    2. INTERPOLATE MEAN VALUES.
!        ------------------------

      IF (EMEAN1.EQ.0.0_JWRB) THEN

!*    2.1 ENERGY OF SPECTRUM 1 IS ZERO.
!         -----------------------------

        EMEAN = GW2*EMEAN2
        FMEAN = FMEAN2
        THETM = THETM2
        DO M=1,ML
          DO K=1,KL
            FL(K,M) = GW2*F2(K,M)
          ENDDO
        ENDDO

        RETURN

      ELSEIF (EMEAN2.EQ.0.0_JWRB) THEN

!*    2.2 ENERGY OF SPECTRUM 2 IS ZERO.
!         -----------------------------

        EMEAN = GW1*EMEAN1
        FMEAN = FMEAN1
        THETM = THETM1
        DO M=1,ML
          DO K=1,KL
            FL(K,M) = GW1*F1(K,M)
          ENDDO
        ENDDO

        RETURN

      ENDIF

!*    2.3 ENERGY BOTH SPECTRUM IS GT ZERO.
!         -------------------------------

      EMEAN = GW1*EMEAN1+GW2*EMEAN2
      FMEAN = GW1*FMEAN1+GW2*FMEAN2
      CM1 = COS(THETM1)
      CM2 = COS(THETM2)
      SM1 = SIN(THETM1)
      SM2 = SIN(THETM2)
      CM  = GW1*CM1 + GW2*CM2
      SM  = GW1*SM1 + GW2*SM2
      THETM = ATAN2(SM,CM)
      THETM = MOD(THETM+ZPI,ZPI)

! ----------------------------------------------------------------------

!*    3. SPECTRUM 1.
!        -----------


!*    3.1 ROTATE.
!         -------


      RTHET1 = THETM - THETM1
      CALL ROTSPEC (NFRE, NANG, ML, KL, F1, F3, RTHET1)

!*    3.2 STRETCH.
!         --------

      GAMMA = FMEAN1/FMEAN
      CALL STRSPEC (NFRE, NANG, ML, KL, FR, F3, GAMMA)

!*    3.3 ADJUST ENERGY.
!         --------------

      EMEANH = EMEAN/EMEAN1
      DO M=1,ML
        DO K=1,KL
          F3(K,M) = F3(K,M)*EMEANH
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    4. SPECTRUM 2.
!        -----------


!*    4.1 ROTATE.
!         -------

      RTHET2 = THETM-THETM2
      CALL ROTSPEC (NFRE, NANG, ML, KL, F2, F4, RTHET2)

!*    4.2 STRETCH.
!         --------

      GAMMA = FMEAN2/FMEAN
      CALL STRSPEC (NFRE, NANG, ML, KL, FR, F4, GAMMA)

!*    4.3 ADJUST ENERGY.
!         --------------

      EMEANH = EMEAN/EMEAN2
      DO M=1,ML
        DO K=1,KL
          F4(K,M) = F4(K,M)*EMEANH
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    5. LINEAR INTERPOLATION OF NEW SPECTRA.
!        ------------------------------------

      DO M=1,ML
        DO K=1,KL
          FL(K,M) = GW1*F3(K,M) + GW2*F4(K,M)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE INTSPEC
