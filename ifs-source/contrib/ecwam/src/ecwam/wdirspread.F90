! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WDIRSPREAD (KIJS, KIJL, F, EMEAN, LLPEAKF, WDIRSPRD)

! ----------------------------------------------------------------------

!**** *WDIRSPREAD* - COMPUTATION OF THE DIRECTIONAL SPREAD 

!     JEAN BIDLOT    ECMWF           MARCH 2000 

!*    PURPOSE.
!     --------

!     TO COMPUTE THE DIRECTIONAL SPREAD EITHER FOR THE WHOLE FREQUENCY
!     RANGE OR BASED ON THE PEAK FREQUENCY.

!**   INTERFACE.
!     ----------

!       *CALL* *WDIRSPREAD (KIJS, KIJL, F, LLPEAKF, WDIRSPRD)*
!          *KIJS*     - INDEX OF FIRST GRIDPOINT.
!          *KIJL*     - INDEX OF LAST GRIDPOINT.
!          *F*        - SPECTRUM.
!          *EMEAN*    - MEAN WAVE ENERGY FOR THE WAVE SYSTEM (INPUT).
!          *LLPEAKF*  - TRUE IF THE DIRECTIONAL SPREAD IS BASED ON 
!                       THE PEAK FREQUENCY.
!          *WDIRSPRD* - DIRECTIONAL SPREAD.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!      *SCOSFL* 
!      *PEAKFRI*

!     REFERENCE.
!     ----------

!       ORIGINALLY REQUESTED BY LUIGI CAVALERI. WE WILL TRY TO PRODUCE
!       SOME FORM OF REFERENCE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,WETAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      USE, INTRINSIC :: IEEE_EXCEPTIONS

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "peakfri.intfb.h"
#include "scosfl.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMEAN
      LOGICAL, INTENT(IN) :: LLPEAKF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: WDIRSPRD


      INTEGER(KIND=JWIM) :: IJ, M
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: IFRINDEX
      REAL(KIND=JWRB) :: COEF_FR, ONE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: F1D 

      LOGICAL :: LL_HALT_INVALID

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDIRSPREAD',0,ZHOOK_HANDLE)

      ! Turn off Floating-Point-Exceptions in this scope to avoid FPE_INVALID in optimized code
      !   with branch prediction. It is safe to do so as DIV_BY_ZERO is protected.
      CALL IEEE_GET_HALTING_MODE(IEEE_INVALID, LL_HALT_INVALID)
      IF (LL_HALT_INVALID)   CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)

!*    1. INITIALIZE ARRAYS
!        -----------------

      ONE=1.0_JWRB
      COEF_FR=WETAIL*FR(NFRE)

      DO IJ = KIJS,KIJL
        WDIRSPRD(IJ) = 0.0_JWRB
      ENDDO

      IF(LLPEAKF) THEN
!     COMPUTATION IS BASED ON THE PEAK FREQUENCY
        CALL PEAKFRI (KIJS, KIJL, F, IFRINDEX, TEMP, F1D)
        CALL SCOSFL (KIJS, KIJL, F, IFRINDEX, WDIRSPRD)
        DO IJ = KIJS,KIJL
          IF(TEMP(IJ) > 0.0_JWRB) THEN
            WDIRSPRD(IJ) = MIN(WDIRSPRD(IJ)/TEMP(IJ),ONE)
          ELSE
            WDIRSPRD(IJ) = ONE 
          ENDIF
        ENDDO

      ELSE
!     COMPUTATION IS BASED ON THE WHOLE FREQUENCY RANGE
        DO M = 1,NFRE
          IFRINDEX=M
          CALL SCOSFL (KIJS, KIJL, F, IFRINDEX, TEMP)
          DO IJ = KIJS,KIJL
            WDIRSPRD(IJ) = WDIRSPRD(IJ) + TEMP(IJ)*DFIM(M)
          ENDDO
        ENDDO
   
!       ADD TAIL CONTRIBUTION
!       note that it uses TEMP because the last call to SCOSFL is made
        DO IJ = KIJS,KIJL
!         the division by delth is needed since dfim contains delth
          WDIRSPRD(IJ) = WDIRSPRD(IJ)/DELTH + TEMP(IJ)*COEF_FR
        ENDDO

!       NORMALISATION
        DO IJ = KIJS,KIJL
          IF(EMEAN(IJ) > EPSMIN) THEN
            WDIRSPRD(IJ) = MIN(WDIRSPRD(IJ)/EMEAN(IJ),ONE)
          ELSE
            WDIRSPRD(IJ) = ONE 
          ENDIF
        ENDDO

      ENDIF

!     COMPUTE THE ACTUAL SPREAD

      DO IJ = KIJS,KIJL
        WDIRSPRD(IJ) = SQRT(2.0_JWRB*(ONE-WDIRSPRD(IJ)))
      ENDDO

      IF (LL_HALT_INVALID)   CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE.)

      IF (LHOOK) CALL DR_HOOK('WDIRSPREAD',1,ZHOOK_HANDLE)

      END SUBROUTINE WDIRSPREAD
