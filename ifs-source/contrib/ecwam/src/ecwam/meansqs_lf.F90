! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MEANSQS_LF(NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS)

! ----------------------------------------------------------------------

!**** *MEANSQS_LF* - COMPUTATION OF RESOLVED MEAN SQUARE SLOPE UP TO WAVE FREQUENCY FR(NFRE_EFF)

!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT FOR THE RESOLVED SPECTRUM UP TO FREQUENCY FR(NFRE_EFF).
!       NFRE_EFF <= NFRE

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS_LF (NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS)*
!              *XKMSS*   - WAVE NUMBER CUT OFF
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *WAVNUM*  - WAVE NUMBER.
!              *XMSS*    - MEAN SQUARE SLOPE (OUTPUT).

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

      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI
      USE YOWFRED  , ONLY : FR       ,ZPIFR    ,DFIM
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE_EFF
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XMSS 


      INTEGER(KIND=JWIM) :: IJ, M, K, KFRE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP1, TEMP2

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MEANSQS_LF',0,ZHOOK_HANDLE)

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      KFRE=MIN(NFRE_EFF, NFRE)

      XMSS(KIJS:KIJL) = 0.0_JWRB

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M = 1, KFRE
          DO IJ = KIJS, KIJL
            TEMP1(IJ) = DFIM(M)*WAVNUM(IJ,M)**2
            TEMP2(IJ) = 0.0_JWRB
          ENDDO
          DO K = 1, NANG
            DO IJ = KIJS, KIJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ = KIJS, KIJL
            XMSS(IJ) = XMSS(IJ)+TEMP1(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS_LF',1,ZHOOK_HANDLE)

      END SUBROUTINE MEANSQS_LF
