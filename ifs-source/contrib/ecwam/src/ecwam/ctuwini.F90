! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CTUWINI (KIJS, KIJL, NINF, NSUP, BLK2GLO, COSPHM1_EXT,   &
 &                  WLATM1, WCORM1, DP)
! ----------------------------------------------------------------------

!**** *CTUWINI* - INITIALISATION FOR CTUW


!*    PURPOSE.
!     --------


! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
USE YOWDRVTYPE  , ONLY : WVGRIDGLO

USE YOWMAP   , ONLY : NGY
USE YOWPARAM , ONLY : NANG     ,NFRE_RED
USE YOWGRID  , ONLY : COSPH
USE YOWTEST  , ONLY : IU06
USE YOWSTAT  , ONLY : ICASE
USE YOWUBUF  , ONLY : KLAT     ,WLAT     ,KCOR     ,WCOR     ,  &
 &                    WLATN    ,WLONN    ,WCORN

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL  ! GRID POINT INDEXES
INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP  ! HALO EXTEND NINF:NSUP+1
TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT ! 1/COSPH
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(OUT) :: WLATM1 ! 1 - WLAT
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,4), INTENT(OUT) :: WCORM1 ! 1 - WCOR
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(OUT) :: DP     ! COS PHI FACTOR


INTEGER(KIND=JWIM) :: IJ, K, M, IC, ICR, ICL, KY, KK, KKM
INTEGER(KIND=JWIM) :: NLAND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTUWINI',0,ZHOOK_HANDLE)

      NLAND = NSUP+1

      DO IC=1,2
        DO IJ = KIJS,KIJL
          IF (KLAT(IJ,IC,1) < NLAND .AND. KLAT(IJ,IC,2) < NLAND) THEN
!           BOTH CLOSEST AND SECOND CLOSEST POINTS ARE OVER THE OCEAN
            WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
          ELSE IF (KLAT(IJ,IC,1) == NLAND) THEN
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           CLOSEST POINT IS OVER LAND
            IF (WLAT(IJ,IC) <= 0.75_JWRB) WLAT(IJ,IC)=0.0_JWRB
            WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
          ELSE
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           SECOND CLOSEST POINT IS OVER LAND
            IF (WLAT(IJ,IC) >= 0.5_JWRB) WLAT(IJ,IC)=1.0_JWRB
            WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
          ENDIF
        ENDDO
      ENDDO

      DO ICR=1,4
        DO IJ = KIJS,KIJL
          IF (KCOR(IJ,ICR,1) < NLAND .AND. KCOR(IJ,ICR,2) < NLAND) THEN
!           BOTH CLOSEST AND SECOND CLOSEST CORNER POINTS ARE OVER THE OCEAN
            WCORM1(IJ,ICR) = 1.0_JWRB - WCOR(IJ,ICR)
          ELSE IF (KCOR(IJ,ICR,1) == NLAND) THEN
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           CLOSEST CORNER POINT IS OVER LAND
            IF (WCOR(IJ,ICR) <= 0.75_JWRB) WCOR(IJ,ICR)=0.0_JWRB
            WCORM1(IJ,ICR) = 1.0_JWRB - WCOR(IJ,ICR)
          ELSE
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           SECOND CLOSEST CORNER POINT IS OVER LAND
            IF (WCOR(IJ,ICR) > 0.5_JWRB) WCOR(IJ,ICR)=1.0_JWRB 
            WCORM1(IJ,ICR) = 1.0_JWRB - WCOR(IJ,ICR)
          ENDIF
        ENDDO
      ENDDO


!     INITIALISATION

      DO ICL=1,2
        DO IC=1,2
          DO M=1,NFRE_RED
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                WLATN(IJ,K,M,IC,ICL)=0.0_JWRB
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO IC=1,2
        DO M=1,NFRE_RED
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              WLONN(IJ,K,M,IC)=0.0_JWRB
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO ICL=1,2
        DO ICR=1,4
          DO M=1,NFRE_RED
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                WCORN(IJ,K,M,ICR,ICL)=0.0_JWRB
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO



      IF (ICASE == 1) THEN

!*      SPHERICAL GRID.
!       ---------------

!*        COMPUTE COS PHI FACTOR FOR ADJOINING GRID POINT.
!         (for all grid points)
          DO IC=1,2
            DO IJ = KIJS,KIJL
              KY=BLK2GLO%KXLT(IJ)
              KK=KY+2*IC-3
              KKM=MAX(1,MIN(KK,NGY))
              DP(IJ,IC) = COSPH(KKM)*COSPHM1_EXT(IJ)
            ENDDO
          ENDDO
       ENDIF

IF (LHOOK) CALL DR_HOOK('CTUWINI',1,ZHOOK_HANDLE)

END SUBROUTINE CTUWINI
