! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MSTART (IOPTI, FETCH, FRMAX, THETAQ,     &
     &                   FM, ALFA, ZGAMMA, SA, SB,        &
     &                   KIJS, KIJL, FL1,                 &
     &                   WSWAVE, WDWAVE) 
! ----------------------------------------------------------------------

!**** *MSTART* - MAKES START FIELDS FOR WAMODEL.

!      H. GUNTHER    ECMWF    MAY 1990
!      H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.
!      J. BIDLOT     ECMWF    FEBRUARY 96  MESSAGE PASSING

!*    PURPOSE.
!     --------

!       TO GENERATE WAMODEL START FIELDS.

!**   INTERFACE.
!     ----------

!   *CALL* *MSTART (IOPTI, FETCH, FRMAX, THETAQ,
!    &              FM, ALFA, ZGAMMA, SA, SB, 
!    &              NPROMA, KIJS, KIJL,
!    &              FL1,WSWAVE,WDWAVE)*
!      *IOPTI*  INTEGER   START FIELD OPTION
!                         = 0 FROM PARAMETERS.
!                         = 1 FROM WINDS CALM ENERGY=0.
!                         = 2 FROM WINDS CALM FROM PARAMETERS.
!      *FETCH*  REAL      FETCH IN METERS.
!      *FRMAX*  REAL      MAXIMUM PEAK FREQUENCY IN HERTZ.
!      *THETAQ* REAL      MEAN DIRECTION (RAD).
!      *ALFA*   REAL      ALPHA PARAMETER
!      *ZGAMMA* REAL      OVERSHOOT FACTOR.
!      *SA*     REAL      LEFT PEAK WIDTH.
!      *SB*     REAL      RIGHT PEAK WIDTH.

!      *FL1*    REAL      2-D SPECTRUM FOR EACH GRID POINT 
!      *WSWAVE* REAL      WIND SPEED. 
!      *WDWAVE* REAL      WIND DIRECTION. 

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *PEAK*      - COMPUTE PARAMETERS FROM WIND FOR A BLOCK.
!       *SPECTRA*   - COMPUTES SPECTRA OF A BLOCK.

!    REFERENCE.
!    ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : DEG
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "peak.intfb.h"
#include "spectra.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IOPTI
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), INTENT(IN) :: FETCH, FRMAX, THETAQ
      REAL(KIND=JWRB), INTENT(IN) :: FM, ALFA, ZGAMMA, SA, SB
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, WDWAVE
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(OUT) :: FL1


      INTEGER(KIND=JWIM) :: M, K, IJ, NGOU
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: FP, ALPHAJ, THES

      CHARACTER(LEN=14), PARAMETER :: ZERO='              '

! ----------------------------------------------------------------------

!*    1. DEFINE SPECTRUM FOR LAND POINTS AND WRITE OUTPUT.
!        -------------------------------------------------

!*    2.1.1 INITIAL VALUES DUE TO OPTION.
!           -----------------------------

        IF (IOPTI == 1) THEN
          DO IJ = KIJS, KIJL
            FP(IJ) = 0.0_JWRB
            ALPHAJ(IJ) = 0.0_JWRB
            THES(IJ) = WDWAVE(IJ)
          ENDDO
        ELSEIF (IOPTI == 0) THEN
          DO IJ = KIJS, KIJL
            FP(IJ) = FM
            ALPHAJ(IJ) = ALFA
            THES(IJ) = THETAQ
          ENDDO
        ELSE
          DO IJ = KIJS, KIJL
            FP(IJ) = FM
            ALPHAJ(IJ) = ALFA
            IF (WSWAVE(IJ) > 0.1E-08_JWRB) THEN
              THES(IJ) = WDWAVE(IJ)
            ELSE
              THES(IJ) = 0.0_JWRB
            ENDIF
          ENDDO
        ENDIF

!*    2.1.2 PEAK FREQUENCY AND ALPHAJONS FROM FETCH LAW.
!           --------------------------------------------
        IF (IOPTI /= 0) THEN
          CALL PEAK (KIJS, KIJL, FETCH, FRMAX, WSWAVE, FP, ALPHAJ)
        ENDIF

!*    2.2 COMPUTE SPECTRA FROM PARAMETERS.
!         --------------------------------
        CALL SPECTRA (KIJS, KIJL, ZGAMMA, SA, SB, FP, ALPHAJ, THES, FL1)

      END SUBROUTINE MSTART
