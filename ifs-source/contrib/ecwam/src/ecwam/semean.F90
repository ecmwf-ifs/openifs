! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SEMEAN (FL1, KIJS, KIJL, EM, LLEPSMIN)

! ----------------------------------------------------------------------

!**** *SEMEAN* - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.

!     S.D. HASSELMANN.
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SEMEAN(FL1, KIJS, KIJL, EM, LLEPSMIN)*
!          *FL1*  - SPECTRUM.
!          *KIJS* - LOCAL INDEX OF FIRST GRIDPOINT
!          *KIJL* - LOCAL  INDEX OF LAST GRIDPOIN
!          *EM*   - MEAN VARIANCE 
!          *LLEPSMIN* - TRUE IF THE WAVE ENERGY IS AT LEAST = EPSMIN

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,WETAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: FL1
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EM
      LOGICAL, INTENT(IN) :: LLEPSMIN


      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB) :: DELT25
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEMEAN',0,ZHOOK_HANDLE)

!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      IF(LLEPSMIN) THEN
        DO IJ=KIJS,KIJL
          EM(IJ) = EPSMIN
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          EM(IJ) = 0.0_JWRB
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=1,NFRE
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
          EM(IJ) = EM(IJ) + DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. ADD TAIL ENERGY.
!        ----------------

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DO IJ=KIJS,KIJL
        EM(IJ) = EM(IJ) + DELT25*TEMP(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SEMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SEMEAN
