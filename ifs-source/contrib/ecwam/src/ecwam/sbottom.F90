! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)

!SHALLOW
! ----------------------------------------------------------------------

!**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.

!     G.J.KOMEN AND Q.D.GAO
!     OPTIMIZED BY L.F. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL

!*    PURPOSE.
!     --------

!       COMPUTATION OF BOTTOM FRICTION DISSIPATION

!**   INTERFACE.
!     ----------

!       *CALL* *SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1*     - SPECTRUM.
!          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*      - TOTAL SOURCE FUNCTION ARRAY
!          *WAVNUM*  - WAVE NUMBER
!          *DEPTH*   - WATER DEPTH

!     METHOD.
!     -------

!       SEE REFERENCES.

!     REFERENCES.
!     -----------

!       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
!       BOUWS AND KOMEN, JPO 13(1983)1653-1658

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE    ,NFRE_RED
      USE YOWPCONS , ONLY : GM1
      USE YOWSHAL  , ONLY : BATHYMAX

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 


      INTEGER(KIND=JWIM):: IJ, K, M
      REAL(KIND=JWRB) :: CONST, ARG
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: SBO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',0,ZHOOK_HANDLE)

      CONST = -2.0_JWRB*0.038_JWRB*GM1
      DO M = 1, NFRE_RED
        DO IJ=KIJS,KIJL
          IF(DEPTH(IJ) < BATHYMAX) THEN
            ARG = 2.0_JWRB* DEPTH(IJ)*WAVNUM(IJ,M)
            ARG = MIN(ARG,50.0_JWRB)
            SBO(IJ,M) = CONST*WAVNUM(IJ,M)/SINH(ARG)
          ELSE
            SBO(IJ,M) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      DO M = 1, NFRE_RED
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ,M)*FL1(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M)+SBO(IJ,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',1,ZHOOK_HANDLE)

      END SUBROUTINE SBOTTOM
