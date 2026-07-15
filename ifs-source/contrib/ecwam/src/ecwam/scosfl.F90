! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SCOSFL (KIJS, KIJL, F, MM, MEANCOSFL)

! ----------------------------------------------------------------------

!**** *SCOSFL* - COMPUTATION OF THE CENTERED FOURIER COEFFICIENT FOR
!                A GIVEN FREQUENCY AT EACH GRID POINT.

!     J.-R. BIDLOT        ECMWF       MARCH 2000

!*    PURPOSE.
!     --------

!     TO COMPUTE THE CENTERED FOURIER COEFFICIENT FOR A GIVEN FREQUENCY 

!**   INTERFACE.
!     ----------

!       *CALL* *SCOSFL (KIJS, KIJL, F, MM, MEANCOSFL)*
!          *KIJS* - INDEX OF FIRST GRIDPOINT
!          *KIJL* - INDEX OF LAST GRIDPOINT
!          *F*    - SPECTRUM.
!          *MM*   - FREQUENCY INDEX (array)
!          *MEANCOSFL* - CENTERED FOURIER COEFFICIENT

!     METHOD.
!     -------

!       INTEGRATION OF SPECTRUM TIMES COS(THETA-MEAN THETA)
!       OVER DIRECTION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWFRED  , ONLY : DELTH    ,TH       ,COSTH    ,SINTH

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: F
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: MEANCOSFL


      INTEGER(KIND=JWIM) :: IJ, K
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: MEANDIR, SI, CI

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SCOSFL',0,ZHOOK_HANDLE)

!*    1. INITIALISE ARRAYS
!        -----------------

      DO IJ=KIJS,KIJL
        SI(IJ) = 0.0_JWRB
        CI(IJ) = 0.0_JWRB
        MEANCOSFL(IJ) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. FIND MEAN DIRECTION
!        -------------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          SI(IJ) = SI(IJ)+SINTH(K)*F(IJ,K,MM(IJ))
          CI(IJ) = CI(IJ)+COSTH(K)*F(IJ,K,MM(IJ))
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        IF (CI(IJ) == 0.0_JWRB .AND. SI(IJ) == 0.0_JWRB) THEN
          MEANDIR(IJ) = 0.0_JWRB
        ELSE
          MEANDIR(IJ) = ATAN2(SI(IJ),CI(IJ))
        ENDIF
      ENDDO

!*    3. COMPUTE FOURRIER COEFFICIENT 
!        -----------------------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          MEANCOSFL(IJ)=MEANCOSFL(IJ)+COS(TH(K)-MEANDIR(IJ))*F(IJ,K,MM(IJ)) 
        ENDDO
      ENDDO
      DO IJ=KIJS,KIJL
        MEANCOSFL(IJ) = DELTH*MEANCOSFL(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SCOSFL',1,ZHOOK_HANDLE)

      END SUBROUTINE SCOSFL
