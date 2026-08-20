! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDISSIP_JAN (KIJS, KIJL, FL1, FLD, SL,   &
     &                        WAVNUM,                     &
     &                        EMEAN, F1MEAN, XKMEAN)

! ----------------------------------------------------------------------

!**** *SDISSIP_JAN* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!                        AUGUST 2020 Added small viscous dissipation term

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP_JAN (KIJS, KIJ, FL1, FLD, SL,
!                            WAVNUM,
!                            EMEAN,F1MEAN, XKMEAN,)*
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *WAVNUM* - WAVE NUMBER
!          *EMEAN*  - MEAN ENERGY DENSITY 
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.

! ---------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DELTH    ,DFIM     ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2
      USE YOWPHYS  , ONLY : CDIS     ,DELTA_SDIS, RNU    ,CDISVIS

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT):: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN):: EMEAN, F1MEAN, XKMEAN


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: SCDFM, CONSD, CONSS, DELTA_SDISM1, CVIS
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: CM, TEMP1, SDS, X
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: XK2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',0,ZHOOK_HANDLE)

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      DELTA_SDISM1=1.0_JWRB-DELTA_SDIS

      CONSS = CDIS*ZPI
      DO IJ=KIJS,KIJL
        SDS(IJ)=CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          X(IJ) = WAVNUM(IJ,M)/XKMEAN(IJ)
          XK2(IJ) = WAVNUM(IJ,M)**2
        ENDDO

        CVIS=RNU*CDISVIS
        DO IJ=KIJS,KIJL
          TEMP1(IJ) = SDS(IJ)*X(IJ)*(DELTA_SDISM1 + DELTA_SDIS*X(IJ)) + CVIS*XK2(IJ)
        ENDDO

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FLD(IJ,K,M) = FLD(IJ,K,M) + TEMP1(IJ)
            SL(IJ,K,M) = SL(IJ,K,M) + TEMP1(IJ)*FL1(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_JAN
