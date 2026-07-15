! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SETWAVPHYS 

! ----------------------------------------------------------------------

! SETS THE WIND INPUT PARAMETERS FOR THE PHYSIC CHOICE (IPHYS)
! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV   ,AFCRV    ,BFCRV,    &
     &                ESH      ,ASH      ,BSH      ,ASWKM    ,BSWKM
USE YOWCOUP  , ONLY : LLGCBZ0  ,LLNORMAGAM
USE YOWPARAM , ONLY : NANG
USE YOWPHYS  , ONLY : BETAMAX  ,ZALP     ,ALPHAMIN ,ALPHA    ,ALPHAPMAX,&
     &                CHNKMIN_U, TAUWSHELTER, TAILFACTOR, TAILFACTOR_PM,&
     &                CDIS     ,DELTA_SDIS, CDISVIS,                    &
     &                DELTA_THETA_RN, RN1_RN, DTHRN_A, DTHRN_U,         &
     &                ANG_GC_A, ANG_GC_B, ANG_GC_C,                     &
     &                SWELLF4,  SWELLF7, SWELLF7M1, Z0TUBMAX, Z0RAT,    &
     &                SSDSC5
USE YOWSTAT  , ONLY : IPHYS
USE YOWTEST  , ONLY : IU06

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETWAVPHYS',0,ZHOOK_HANDLE)

      IF (IPHYS == 0) THEN
!       JANSSSEN WIND INPUT PHYSICS:
        ZALP    = 0.008_JWRB
        TAILFACTOR = 2.5_JWRB
        ALPHAMIN = 0.0001_JWRB
        ALPHAPMAX = 0.03_JWRB
        TAUWSHELTER = 0.0_JWRB

!       ANGULAR ADJUSTMENT PARAMETERS FOR THE GRAVITY-CAPILLARY MODEL
        IF (NANG <= 24) THEN
          ANG_GC_A = 0.40_JWRB
          ANG_GC_B = 0.60_JWRB
          ANG_GC_C = 3.0_JWRB
        ELSE
          ANG_GC_A = 0.35_JWRB
          ANG_GC_B = 0.65_JWRB
          ANG_GC_C = 3.0_JWRB
        ENDIF

!       DIRECTIONALITY CORRECTION FACTORS IN THE GOWTH RATE RENORMALISATION 
        DELTA_THETA_RN = 0.75_JWRB
        DTHRN_A = 0.80_JWRB
        DTHRN_U = 33.0_JWRB

!       DIRECTIONALITY CORRECTION FACTOR FOR THE GRAVITY-CAPILLARY MODEL
        RN1_RN = 0.25_JWRB

        TAILFACTOR_PM = 0.0_JWRB   ! i.e. not used


        IF(LLGCBZ0) THEN
          !!! not yet fully tested !!!
          ALPHA   = 0.0055_JWRB
          CHNKMIN_U = 28._JWRB

          IF(LLNORMAGAM) THEN
            BETAMAX = 1.32_JWRB
          ELSE
            !!! untested
            BETAMAX = 1.25_JWRB
          ENDIF

          CDIS = -1.3_JWRB
          DELTA_SDIS = 0.6_JWRB
          CDISVIS = -4.0_JWRB

        ELSE
          ALPHA   = 0.0065_JWRB
          CHNKMIN_U = 33._JWRB
          BETAMAX = 1.20_JWRB

          CDIS = -1.33_JWRB
          DELTA_SDIS = 0.5_JWRB
          CDISVIS = 0.0_JWRB
        ENDIF

!!!     EMPIRICAL CONSTANCE FOR  SPECTRAL UPDATE FOLLOWING DATA ASSIMILATION
        EGRCRV = 1108.0_JWRB
        AGRCRV = 0.06E+6_JWRB
        BGRCRV = 9.70_JWRB
        AFCRV = 4.0E-4_JWRB
        BFCRV = -3.0_JWRB
        ESH = 1711.0_JWRB
        ASH = 8.0E-4_JWRB 
        BSH = 0.96_JWRB 
        ASWKM = 0.0981_JWRB
        BSWKM = 0.425_JWRB


      ELSE IF (IPHYS == 1) THEN
!       ARDHUIN ET AL. (2010) WIND INPUT PHYSICS
        ZALP    = 0.008_JWRB
        TAILFACTOR = 2.5_JWRB
        TAILFACTOR_PM = 3.0_JWRB
        ALPHAMIN = 0.0001_JWRB

!       ANGULAR ADJUSTMENT PARAMETERS FOR THE GRAVITY-CAPILLARY MODEL
        IF (NANG <= 24) THEN
          ANG_GC_A = 0.40_JWRB
          ANG_GC_B = 0.60_JWRB
          ANG_GC_C = 3.0_JWRB
        ELSE
          ANG_GC_A = 0.35_JWRB
          ANG_GC_B = 0.65_JWRB
          ANG_GC_C = 3.0_JWRB
        ENDIF

!       directionality correction factors in the gowth rate renormalisation 
        DELTA_THETA_RN = 0.75_JWRB
        DTHRN_A = 0.60_JWRB
        DTHRN_U = 33.0_JWRB

!       DIRECTIONALITY CORRECTION FACTOR FOR THE GRAVITY-CAPILLARY MODEL
        RN1_RN = 0.25_JWRB

        IF(LLGCBZ0) THEN
          ALPHA   = 0.0055_JWRB
          CHNKMIN_U = 28._JWRB

          SWELLF4 = 1.15E05_JWRB
          SWELLF7 = 4.32E05_JWRB
          SWELLF7M1 = 1.0_JWRB/SWELLF7 

          SSDSC5  = 0.0_JWRB

          Z0TUBMAX = 0.05_JWRB
          Z0RAT = 0.02_JWRB

          ALPHAPMAX = 0.03_JWRB

          IF(LLNORMAGAM) THEN
            BETAMAX = 1.39_JWRB
            TAUWSHELTER = 0.0_JWRB
          ELSE
           !!! not yet fully tested !!!
            BETAMAX = 1.44_JWRB
            TAUWSHELTER = 0.25_JWRB
          ENDIF

        ELSE 
          ALPHA   = 0.0065_JWRB
          CHNKMIN_U = 33._JWRB

          SWELLF4 = 1.5E05_JWRB
          SWELLF7 = 3.6E05_JWRB
          SWELLF7M1 = 1.0_JWRB/SWELLF7 

          SSDSC5  = 0.0_JWRB

          Z0TUBMAX = 0.0005_JWRB
          Z0RAT = 0.04_JWRB

          ALPHAPMAX = 0.031_JWRB

          IF(LLNORMAGAM) THEN
            !!! not yet tested !!!
            BETAMAX = 1.40_JWRB
            TAUWSHELTER = 0.0_JWRB
          ELSE
            BETAMAX = 1.40_JWRB
            TAUWSHELTER = 0.25_JWRB
          ENDIF
        ENDIF

!!!     EMPIRICAL CONSTANCE FOR  SPECTRAL UPDATE FOLLOWING DATA ASSIMILATION
        EGRCRV = 1065.0_JWRB
        AGRCRV = 0.0655E+6_JWRB
        BGRCRV =  10.906_JWRB
        AFCRV = 2.453E-4_JWRB
        BFCRV = -3.1236_JWRB
        ESH = 1711.0_JWRB
        ASH = 8.0E-4_JWRB 
        BSH = 0.96_JWRB 
        ASWKM = 0.0981_JWRB
        BSWKM = 0.425_JWRB

      ELSE
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR IN SETWAVPHYS                 *'
        WRITE (IU06,*) '*  UKNOWN PHYSICS SELECTION :       *'
        WRITE (IU06,*) '*  IPHYS =' , IPHYS
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

IF (LHOOK) CALL DR_HOOK('SETWAVPHYS',1,ZHOOK_HANDLE)

END SUBROUTINE SETWAVPHYS
