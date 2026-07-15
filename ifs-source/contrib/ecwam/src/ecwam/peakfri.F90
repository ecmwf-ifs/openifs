! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE PEAKFRI (KIJS, KIJL, F, IPEAKF, EPEAKF, F1D)

! ----------------------------------------------------------------------

!**** *PEAKFRI* - COMPUTATION OF INDEX OF THE PEAK FREQUENCY AND
!                 THE CORRESPONDING VALUE OF THE 1-D SPECTRAM
!                 AT EACH GRID POINT.

!     JEAN BIDLOT    ECMWF           MARCH 2000 

!*    PURPOSE.
!     --------

!       COMPUTE PEAK FREQUENCY INDEX AT EACH GRID POINT.
!       IF NOT PEAK FOUND THAN EPEAKF WILL BE SET TO 0. AND
!       IPEAKF TO NFRE

!**   INTERFACE.
!     ----------

!       *CALL* *PEAKFRI (KIJS, KIJL, F, IPEAKF, EPEAKF, F1D)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *F*       - SPECTRUM.
!          *IPEAKF   - INDEX OF PEAK FREQUENCY 
!          *EPEAKF*  - 1-D SPECTRUM AT PEAK FREQUENCY
!          *F1D*     - 1D FREQUENCY SPECTRUM

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : DELTH    ,FR
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN) :: F(KIJS:KIJL,NANG,NFRE)
      INTEGER(KIND=JWIM), INTENT(OUT) :: IPEAKF(KIJS:KIJL) 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EPEAKF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: F1D

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PEAKFRI',0,ZHOOK_HANDLE)

!*    1. INITIALIZE ARRAYS
!        -----------------

      DO IJ = KIJS,KIJL
        EPEAKF(IJ) = 0.0_JWRB
        IPEAKF(IJ) = NFRE 
      ENDDO

!*    2. LOOP OVER FREQUENCIES
!        ---------------------

      DO M = 1, NFRE

!*    2.1 COMPUTE 1-D SPECTRUM
!        ---------------------

        DO IJ = KIJS,KIJL
          F1D(IJ,M) = 0.0_JWRB
        ENDDO
        DO K = 1,NANG
          DO IJ = KIJS,KIJL
            F1D(IJ,M) = F1D(IJ,M) + F(IJ,K,M)*DELTH
          ENDDO
        ENDDO

!*    2.2 DEFINE PEAK FREQUENCY
!         ---------------------

        DO IJ = KIJS,KIJL
          IF (EPEAKF(IJ) < F1D(IJ,M)) THEN
            EPEAKF(IJ) = F1D(IJ,M)
            IPEAKF(IJ) = M 
          ENDIF
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('PEAKFRI',1,ZHOOK_HANDLE)

      END SUBROUTINE PEAKFRI
