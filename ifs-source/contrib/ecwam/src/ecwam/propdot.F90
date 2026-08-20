! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE PROPDOT(KIJS, KIJL, NINF, NSUP,                &
     &                   BLK2GLO,                               &
     &                   WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
     &                   COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT,  & 
     &                   THDC, THDD, SDOT)

! ----------------------------------------------------------------------

!**** *PROPDOT* - PROPAGATION DOT TERMS FROM DEPTH AND CURRENT GRADIENT.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF        FEBRUARY 97 MESSAGE PASSING
!     J. BIDLOT    ECMWF        JANUARY 2004 REMOVE MULTI-BLOCK OPTION.

!*    PURPOSE.
!     --------

!       COMPUTATION OF COMMON REFDOT FOR PROPAGATION.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPDOT*(KIJS, KIJL, NINF, NSUP,
!                        BLK2GLO,
!                        WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
!                        DEPTH_EXT, U_EXT, V_EXT,       & 
!                        THDC, THDD, SDOT)
!          *KIJS*        - STARTING INDEX
!          *KIJL*        - ENDING INDEX
!          *NINF:NSUP+1* - 1st DIMENSION OF *_EXT ARRAYS 
!          *BLK2GLO*     - BLOCK TO GRID TRANSFORMATION
!          *WAVNUM_EXT*  - WAVE NUMBER.
!          *CGROUP_EXT*  - GROUP SPPED.
!          *OMOSNH2KD_EXT*- OMEGA / SINH(2KD)
!          *DEPTH_EXT*    - WATER DEPTH (including the halo points)
!          *U_EXT*        - U-COMPONENT OF SURFACE CURRENT (including the halo points)
!          *V_EXT*        - V-COMPONENT OF SURFACE CURRENT (including the halo points)
!         THDC and THDD ARRAYS TO KEEP DEPTH AND CURRENT REFRACTION FOR THETA DOT
!         AND SDOT FOR SIGMA DOT

!     METHOD.
!     -------

!       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED

!     EXTERNALS.
!     ----------

!       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWMAP   , ONLY : NIBLO
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED
      USE YOWSTAT  , ONLY : ICASE    ,IREFRA

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "gradi.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, NINF, NSUP
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(INOUT) :: THDC, THDD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE_RED), INTENT(INOUT) :: SDOT


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: CD, SD, SS, SC, CC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DDPHI, DDLAM, DUPHI, DULAM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DVPHI, DVLAM, DCO, OMDD

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PROPDOT',0,ZHOOK_HANDLE)


!*    2.2 DEPTH AND CURRENT GRADIENTS.
!         ----------------------------

        CALL GRADI (KIJS, KIJL, NINF, NSUP, IREFRA, &
     &              BLK2GLO,                        &
     &              DEPTH_EXT, U_EXT, V_EXT,        & 
     &              DDPHI, DDLAM, DUPHI,            &
     &              DULAM, DVPHI, DVLAM)


!*    2.3 COSINE OF LATITUDES IF SPHERICAL PROPAGATION.
!         ---------------------------------------------

        IF (ICASE == 1) THEN
          DO IJ = KIJS,KIJL
            DCO(IJ) = COSPHM1_EXT(IJ)
          ENDDO
        ELSE
          DO IJ = KIJS,KIJL
            DCO(IJ) = 1.0_JWRB
          ENDDO
        ENDIF

!*    2.4 DEPTH GRADIENT PART OF SIGMA DOT.
!         ---------------------------------

        IF (IREFRA == 3) THEN
          DO IJ = KIJS,KIJL
            OMDD(IJ) = V_EXT(IJ)*DDPHI(IJ) + U_EXT(IJ)*DDLAM(IJ)*DCO(IJ)
          ENDDO
        ELSEIF (IREFRA == 2) THEN
          DO IJ = KIJS,KIJL
            OMDD(IJ) = 0.0_JWRB
          ENDDO
        ENDIF


!*    2.5. LOOP OVER DIRECTIONS.
!          ---------------------

        DO K=1,NANG
          SD = SINTH(K)
          CD = COSTH(K)

!*    2.5.1. DEPTH GRADIENT OF THETA DOT.
!            ----------------------------

          IF (IREFRA == 1 .OR. IREFRA == 3) THEN
            DO IJ = KIJS,KIJL
              THDD(IJ,K) = SD*DDPHI(IJ) - CD*DDLAM(IJ)*DCO(IJ)
            ENDDO
          ELSE
            DO IJ = KIJS,KIJL
              THDD(IJ,K) = 0.0_JWRB
            ENDDO
          ENDIF

!*    2.5.2 SIGMA DOT AND THETA DOT PART FROM CURRENT GRADIENT.
!           ---------------------------------------------------

          IF (IREFRA == 2 .OR. IREFRA == 3) THEN

            SS  = SD**2
            SC  = SD*CD
            CC  = CD**2
            DO IJ = KIJS,KIJL
              SDOT(IJ,K,NFRE_RED) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)      &
     &                        - (SS*DULAM(IJ) + SC*DVLAM(IJ))*DCO(IJ)
              THDC(IJ,K) =  SS*DUPHI(IJ) + SC*DVPHI(IJ)               &
     &                    - (SC*DULAM(IJ) + CC*DVLAM(IJ))*DCO(IJ)
            ENDDO

!*    2.5.3 LOOP OVER FREQUENCIES.
!           ----------------------

            DO M=1,NFRE_RED
              DO IJ=KIJS,KIJL
                SDOT(IJ,K,M) = (SDOT(IJ,K,NFRE_RED)*CGROUP_EXT(IJ,M)   &
     &           + OMDD(IJ)*OMOSNH2KD_EXT(IJ,M)) * WAVNUM_EXT(IJ,M) 
              ENDDO
            ENDDO

          ENDIF

!*      BRANCH BACK TO 2.5 FOR NEXT DIRECTION.
        ENDDO

      IF (LHOOK) CALL DR_HOOK('PROPDOT',1,ZHOOK_HANDLE)

      END SUBROUTINE PROPDOT
