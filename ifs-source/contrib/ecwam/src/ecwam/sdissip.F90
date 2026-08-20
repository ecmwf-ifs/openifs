! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDISSIP (KIJS, KIJL, FL1, FLD, SL,  &
     &                    WAVNUM, CGROUP, XK2CG,     &
     &                    EMEAN, F1MEAN, XKMEAN,     &
     &                    UFRIC, COSWDIF, RAORW)
! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (KIJS, KIJL, FL1, FLD, SL, *
!                        WAVNUM, CGROUP, XK2CG,  
!                        EMEAN, F1MEAN, XKMEAN,*
!                        UFRIC, COSWDIF, RAORW)*
!         *KIJS* - INDEX OF FIRST GRIDPOINT
!         *KIJL* - INDEX OF LAST GRIDPOINT
!         *FL1*  - SPECTRUM.
!         *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!       *WAVNUM* - WAVE NUMBER
!       *CGROUP* - GROUP SPEED
!       *XK2CG*  - (WAVNUM)**2 * GROUP SPEED
!       *EMEAN* - MEAN ENERGY DENSITY 
!       *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!       *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!       *UFRIC*  - FRICTION VELOCITY IN M/S.
!       *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
!       *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "sdissip_ard.intfb.h"
#include "sdissip_jan.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, RAORW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0)
         CALL SDISSIP_JAN (KIJS, KIJL, FL1 ,FLD, SL,  &
     &                     WAVNUM,                    &
     &                     EMEAN, F1MEAN, XKMEAN)

      CASE(1) 
         CALL SDISSIP_ARD (KIJS, KIJL, FL1 ,FLD, SL,   &
     &                     WAVNUM, CGROUP, XK2CG,      &
     &                     UFRIC, COSWDIF, RAORW)
      END SELECT 

      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP
