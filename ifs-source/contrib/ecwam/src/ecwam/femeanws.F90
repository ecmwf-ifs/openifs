! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE FEMEANWS (KIJS, KIJL, FL1, XLLWS, EM, FM)

! ----------------------------------------------------------------------

!**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY 
!                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
!                  BY XLLWS

!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
!       SPECTRUM WHERE XLLWS IS NON ZERO.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEANWS (KIJS, KIJL, FL1, XLLWS, EM, FM)*
!              *KIJS*   - INDEX OF FIRST GRIDPOINT
!              *KIJL*   - INDEX OF LAST GRIDPOINT
!              *FL1*    - SPECTRUM.
!              *XLLWS* - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
!              *EM*     - MEAN WAVE ENERGY (OUTPUT)
!              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,    &
     &                WETAIL    ,FRTAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(kIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1, XLLWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EM, FM


      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB) :: DELT25, DELT2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FEMEANWS',0,ZHOOK_HANDLE)

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=KIJS,KIJL
        EM(IJ) = EPSMIN
        FM(IJ) = EPSMIN
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH


!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------
      
      DO M=1,NFRE
        TEMP2(KIJS:KIJL) = 0.0_JWRB
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            TEMP2(IJ) = TEMP2(IJ)+XLLWS(IJ,K,M)*FL1(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EM(IJ) = EM(IJ)+DFIM(M)*TEMP2(IJ)
          FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
        ENDDO
      ENDDO

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      DO IJ=KIJS,KIJL
        EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
        FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
        FM(IJ) = EM(IJ)/FM(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FEMEANWS',1,ZHOOK_HANDLE)

      END SUBROUTINE FEMEANWS
