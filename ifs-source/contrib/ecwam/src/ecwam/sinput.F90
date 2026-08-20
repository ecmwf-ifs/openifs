! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SINPUT (NGST, LLSNEG, KIJS, KIJL, FL1, & 
     &                   WAVNUM, CINV, XK2CG,           &
     &                   WDWAVE, WSWAVE, UFRIC, Z0M,    &
     &                   COSWDIF, SINWDIF2,             & 
     &                   RAORW, WSTAR, RNFAC,           &
     &                   FLD, SL, SPOS, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.


!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (NGST, LLSNEG, KIJS, KIJL, FL1,
!    &                WAVNUM, CINV, XK2CG,
!    &                WDWAVE, UFRIC, Z0M,
!    &                COSWDIF, SINWDIF2,
!    &                RAORW, WSTAR, FLD, SL, SPOS, XLLWS)
!         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
!                - IF = 2 THEN GUSTINESS PARAMETERISATION
!         *LLSNEG* - IF TRUE THEN THE NEGATIVE SINPUT WILL BE COMPUTED
!         *KIJS* - INDEX OF FIRST GRIDPOINT.
!         *KIJL* - INDEX OF LAST GRIDPOINT.
!          *FL1* - SPECTRUM.
!       *WAVNUM* - WAVE NUMBER.
!         *CINV* - INVERSE PHASE VELOCITY.
!       *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPPED.
!       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *UFRIC* - FRICTION VELOCITY IN M/S.
!        *Z0M*   - ROUGHNESS LENGTH IN M.
!      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
!     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
!        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
!        *XLLWS* - = 1 WHERE SINPUT IS POSITIVE

!     METHOD.
!     -------

!     DEPENDING ON THE VALUE OF IPHYS, DIFFERENT INPUTE SOURCE TERm WILL BE CALLED

!     EXTERNALS.
!     ----------

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "sinput_ard.intfb.h"
#include "sinput_jan.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      LOGICAL, INTENT(IN) :: LLSNEG
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, XK2CG

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, WSWAVE, UFRIC, Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: RAORW, WSTAR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF, SINWDIF2

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0)
        CALL SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1,  &
     &                   WAVNUM, CINV, XK2CG,            &
     &                   WSWAVE, UFRIC, Z0M,     &
     &                   COSWDIF, SINWDIF2,              & 
     &                   RAORW, WSTAR, RNFAC,            &
     &                   FLD, SL, SPOS, XLLWS)
      CASE(1) 
        CALL SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1,  &
     &                   WAVNUM, CINV, XK2CG,            &
     &                   WDWAVE, WSWAVE, UFRIC, Z0M,     &
     &                   COSWDIF, SINWDIF2,              & 
     &                   RAORW, WSTAR, RNFAC,            &
     &                   FLD, SL, SPOS, XLLWS)
      END SELECT 

      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT
