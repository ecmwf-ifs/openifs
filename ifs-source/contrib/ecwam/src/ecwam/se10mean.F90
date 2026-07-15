! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SE10MEAN (KIJS, KIJL, FL1, E10)

! ----------------------------------------------------------------------

!**** *SE10MEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                 WITH PERIOD LARGER THAN 10s

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SE10MEAN(KIJS, KIJL, FL1, E10)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1      - SPECTRUM.
!          *E10*     - MEAN ENERGY

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: E10


      INTEGER(KIND=JWIM) :: IJ, M, K, MCUT

      REAL(KIND=JWRB), PARAMETER :: FCUT=0.1_JWRB
      REAL(KIND=JWRB) :: DFCUT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFIMLOC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',0,ZHOOK_HANDLE)

      MCUT=NINT(LOG10(FCUT/FR(1))*FLOGSPRDM1)+1
      MCUT=MIN(MAX(1,MCUT),NFRE)

      DFIMLOC(:)=DFIM(:)
      DFCUT=0.5_JWRB*MAX(0.0_JWRB,FCUT-0.5_JWRB*(FR(MAX(1,MCUT-1))+FR(MCUT)))
      DFIMLOC(MCUT)=MIN(DFCUT*DELTH,DFIM(MCUT))


!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=KIJS,KIJL
        E10(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=1,MCUT
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          E10(IJ) = E10(IJ)+DFIMLOC(M)*TEMP(IJ)
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SE10MEAN
