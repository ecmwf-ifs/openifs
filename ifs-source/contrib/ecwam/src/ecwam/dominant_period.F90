! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE DOMINANT_PERIOD (KIJS, KIJL, FL1, DP)

! ----------------------------------------------------------------------

!**** *DOMINANT_PERIOD* - COMPUTATION OF THE DOMINANT PERIOD

!*    PURPOSE.
!     --------

!       COMPUTE DOMINANT PERIOD AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *DOMINANT_PERIOD (KIJS, KIJL, FL1, DP)*
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *FL1*     - 2D-SPECTRUM.
!              *DP*      - DOMINANT PERIOD

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMFR   ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK   ,JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: DP

      REAL(KIND=JWRB), PARAMETER :: FLTHRS = 0.1_JWRB

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP, EM, FCROP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: F1D4

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('DOMINANT_PERIOD',0,ZHOOK_HANDLE)

      F1D4(:,:) = 0.0_JWRB
      EM(:) = 0.0_JWRB
      DP(:) = 0.0_JWRB
      FCROP(:) = 0.0_JWRB

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            IF( FL1(IJ,K,M) > FCROP(IJ) ) THEN
               FCROP(IJ) = FL1(IJ,K,M)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      FCROP(:) = FLTHRS*FCROP(:)

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            IF(FL1(IJ,K,M) > FCROP(IJ)) F1D4(IJ,M) = F1D4(IJ,M)+FL1(IJ,K,M)*DELTH
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          F1D4(IJ,M) = F1D4(IJ,M)**4
          EM(IJ) = EM(IJ)+DFIM(M)*F1D4(IJ,M)
          DP(IJ) = DP(IJ)+DFIMFR(M)*F1D4(IJ,M)
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        IF(EM(IJ).GT.0.0_JWRB .AND. DP(IJ).GT.EPSMIN ) THEN
          DP(IJ) = EM(IJ)/DP(IJ)
        ELSE
          DP(IJ) = 0.0_JWRB 
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('DOMINANT_PERIOD',1,ZHOOK_HANDLE)

      END SUBROUTINE DOMINANT_PERIOD
