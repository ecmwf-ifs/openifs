! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MWP1 (KIJS, KIJL, F, MEANWP1)

! ----------------------------------------------------------------------

!**** *MWP1* - COMPUTATION OF MEAN PERIOD BASED ON THE FIRST MOMENT
!              FOR EACH GRID POINT.

!     J-R BIDLOT    ECMWF     MARCH 2000

!*    PURPOSE.
!     --------

!       COMPUTE MEAN PERIOD 1 AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MWP1 (KIJS, KIJL, F, MEANWP1)*
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *MEANWP1* - MEAN PERIOD BASED ON THE FIRST MOMENT.

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

      USE YOWFRED  , ONLY : FR       ,DFIM_SIM ,DFIMFR_SIM,DELTH  ,WETAIL, WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      USE, INTRINSIC :: IEEE_EXCEPTIONS

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN) :: F(KIJS:KIJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: MEANWP1


      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELT25, COEF_FR, FR1M1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP, EM

      LOGICAL :: LL_HALT_INVALID

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MWP1',0,ZHOOK_HANDLE)

      ! Turn off Floating-Point-Exceptions in this scope to avoid FPE_INVALID in optimized code
      !   with branch prediction. It is safe to do so as DIV_BY_ZERO is protected.
      CALL IEEE_GET_HALTING_MODE(IEEE_INVALID, LL_HALT_INVALID)
      IF (LL_HALT_INVALID)   CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)

      DO IJ=KIJS,KIJL
        EM(IJ) = 0.0_JWRB
        MEANWP1(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE_ODD
        DO IJ=KIJS,KIJL
          TEMP(IJ) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EM(IJ) = EM(IJ)+DFIM_SIM(M)*TEMP(IJ)
          MEANWP1(IJ) = MEANWP1(IJ)+DFIMFR_SIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      FR1M1 = 1.0_JWRB/FR(1)
      DELT25 = WETAIL*FR(NFRE_ODD)*DELTH
      COEF_FR = WP1TAIL*DELTH*FR(NFRE_ODD)**2

      DO IJ=KIJS,KIJL
        EM(IJ) = EM(IJ)+DELT25*TEMP(IJ)
        MEANWP1(IJ) = MEANWP1(IJ)+COEF_FR*TEMP(IJ)
        IF(EM(IJ) > 0.0_JWRB .AND. MEANWP1(IJ) > EPSMIN ) THEN
          MEANWP1(IJ) = EM(IJ)/MEANWP1(IJ)
          MEANWP1(IJ) = MIN(MEANWP1(IJ),FR1M1)
        ELSE
          MEANWP1(IJ) = 0.0_JWRB
        ENDIF
      ENDDO

      IF (LL_HALT_INVALID)   CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE.)

      IF (LHOOK) CALL DR_HOOK('MWP1',1,ZHOOK_HANDLE)

      END SUBROUTINE MWP1
