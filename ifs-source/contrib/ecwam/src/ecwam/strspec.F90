! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE STRSPEC (NFRE, NANG, ML, KL, FR, FL, GAMMA)

! ----------------------------------------------------------------------

!**** *STRSPEC* - ROUTINE TO STRETCH A SPECTRUM.

!      EVA BAUER      MPI  HAMBURG    MAY 1990.
!      H. GUNTHER     GKSS/ECMWF      JAN 1991  MODIFIED FOR CYCLE_4.

!*     PURPOSE.
!      -------


!**   INTERFACE.
!     ----------

!       *CALL* *STRSPEC (NFRE, NANG, ML, KL, FR, FL, GAMMA)*
!         *NFRE*  - FREQUENCY DIMENSION OF ARRAYS.
!         *NANG*  - DIRECTION DIMENSION OF ARRAYS.
!         *ML*    - NUMBER OF FREQUENCIES.
!         *KL*    - NUMBER OF DIRECTIONS.
!         *FR*    - FREQUENCY ARRAY.
!         *FL*    - INPUT AND OUTPUT SPECTRUM.
!         *GAMMA* - STRETCHING PARAMETER.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ---------

!       NONE.

!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG, ML, KL
      REAL(KIND=JWRB), INTENT(IN) :: GAMMA
      REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(IN) :: FR
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(INOUT) :: FL


      INTEGER(KIND=JWIM) :: M, K
      INTEGER(KIND=JWIM) :: INC, MC, IFR, IFRP1
      REAL(KIND=JWRB) :: ALO, GAMS, Z, ADIF, BDIF
      REAL(KIND=JWRB), DIMENSION(ML) :: AR2
      REAL(KIND=JWRB), DIMENSION(KL,ML) :: AR1

!---------------------------------------------------------------------

      IF (GAMMA.EQ.1.0_JWRB) RETURN

!*    1. INITIALIZATION.
!        ---------------

      DO M=1,ML
        DO K=1,KL
          AR1(K,M) = 0.0_JWRB
        ENDDO
      ENDDO
!!!!! 1.1 should actually be FRATIO !!!!
      ALO = LOG10(1.1_JWRB)
      GAMS = GAMMA

!*    2. DETERMINE ACROSS HOW MANY FREQUENCY BINS THE
!        STRETCHING IS ACTING AND THE STRECHED FREQUENCIES.
!        ---------------------------------------------------

      INC = INT(LOG10(GAMS)/ALO)
!!!!! 1.1 should actually be FRATIO !!!!
      Z = ABS(1.1_JWRB**INC - GAMS)
      DO M=1,ML
        AR2(M) = FR(M) * GAMS
      ENDDO

!*    3. STRECH SPECTRUM.
!        ----------------

      IF (Z.LE.0.001_JWRB) THEN

!*    3.1 SHIFT SPECTRUM IF GAMMA IS A POWER OF 1.1.
!         ------------------------------------------

        IF (GAMS.GT.1.0_JWRB) THEN

!*    3.1.1 SHIFT TO LOWER FREQUENIES.
!           --------------------------

          DO M = 1,ML-INC
            MC = M + INC
            DO K = 1,KL
              AR1(K,M) = FL(K,MC)
            ENDDO
          ENDDO
        ELSE

!*    3.1.2 SHIFT TO HIGHER FREQUENCIES.
!           ----------------------------

          DO M = 1-INC,ML
            MC = M + INC
            DO K = 1,NANG
              AR1(K,M) = FL(K,MC)
            ENDDO
          ENDDO
        ENDIF
      ELSE

!*    3.2 SHIFT AND LINEAR INTERPOLATION OF SPECTRAL ENERGY
!*        IF GAMMA IS NOT A POWER OF 1.1.
!*
!*        (SPECTRUM HAS ZERO ENERGY AT FREQUENCY FR(NFRE) )
!         -------------------------------------------------

        IF (GAMS.GT.1.0_JWRB) THEN

!*    3.2.1 SHIFT TO LOWER FREQUENCIES.
!           ---------------------------

          DO M=1,ML-INC-1
            IFR = INT(LOG10(AR2(M)/FR(1))/ALO+1.0_JWRB)
            IFRP1 = IFR+1
            MC = M + INC
            ADIF = (FR(IFRP1)-AR2(M)) / (FR(IFRP1)-FR(IFR))
            BDIF = 1. - ADIF
            DO K=1,KL
              AR1(K,M) = ADIF * FL(K,MC) + BDIF * FL(K,MC+1)
            ENDDO
          ENDDO
        ELSE

!*     3.2.2 SHIFT TO HIGHER FREQUENCIES.
!            ----------------------------

          DO M = 2-INC,ML
            IFR = INT(LOG10(AR2(M)/FR(1))/ALO+1.0_JWRB)
            IFRP1 = IFR+1
            MC = M + INC - 1
            ADIF = (FR(IFRP1)-AR2(M)) / (FR(IFRP1)-FR(IFR))
            BDIF = 1. - ADIF
            DO K=1,KL
              AR1(K,M) = ADIF * FL(K,MC) + BDIF * FL(K,MC+1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!*    4. COPY NEW TO OLD.
!        ----------------


      DO M=1,ML
        DO K=1,KL
          FL(K,M) = AR1(K,M)
        ENDDO
      ENDDO

      END SUBROUTINE STRSPEC
