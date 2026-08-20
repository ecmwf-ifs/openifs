! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE STRESSO (KIJS, KIJL, MIJ, RHOWGDFTH,         &
     &                    FL1, SL, SPOS,                      &
     &                    CINV,                               &
     &                    WDWAVE, UFRIC, Z0M, AIRD, RNFAC,    &
     &                    COSWDIF, SINWDIF2,                  &
     &                    TAUW, TAUWDIR, PHIWA, LLPHIWA)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
!     J. BIDLOT            ECMWF FEBRUARY   1996-97
!     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!        *CALL* *STRESSO (KIJS, KIJL, MIJ, RHOWGDFTH, 
!                         FL1, SL, SPOS,
!    &                    CINV,
!    &                    WDWAVE, UFRIC, Z0M, AIRD, RNFAC,
!    &                    COSWDIF, SINWDIF2,
!    &                    TAUW, TAUWDIR, PHIWA)*
!         *KIJS*        - INDEX OF FIRST GRIDPOINT.
!         *KIJL*        - INDEX OF LAST GRIDPOINT.
!         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
!         *FL1*         - WAVE SPECTRUM.
!         *SL*          - WIND INPUT SOURCE FUNCTION ARRAY (positive and negative contributions).
!         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
!         *CINV*        - INVERSE PHASE VELOCITY.
!         *WDWAVE*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                         NOTATION (POINTING ANGLE OF WIND VECTOR,
!                         CLOCKWISE FROM NORTH).
!         *UFRIC*       - FRICTION VELOCITY IN M/S.
!         *Z0M*         - ROUGHNESS LENGTH IN M.
!         *AIRD*        - AIR DENSITY IN KG/M**3.
!         *RNFAC*       - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!         *COSWDIF*     - COS(TH(K)-WDWAVE(IJ))
!         *SINWDIF2*    - SIN(TH(K)-WDWAVE(IJ))**2
!         *TAUW*        - KINEMATIC WAVE STRESS IN (M/S)**2
!         *TAUWDIR*     - KINEMATIC WAVE STRESS DIRECTION
!         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                         OVER THE FULL FREQUENCY RANGE.
!         *LLPHIWA*     - TRUE IF PHIWA NEEDS TO BE COMPUTED

!     METHOD.
!     -------

!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
!       AND DIRECTIONS.
!       BECAUSE ARRAY *SPOS* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SPOS* (CALL FIRST SINPUT, THEN
!       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)

!     REFERENCE.
!     ----------
!       P. JANSSEN, 

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWFRED  , ONLY : FR       ,RHOWG_DFIM ,DELTH    ,TH       ,    &
     &            COSTH    ,SINTH    ,FR5
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPHYS  , ONLY : TAUWSHELTER
      USE YOWTABL  , ONLY : EPS1
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "tau_phi_hf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJS:KIJL)
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, UFRIC, Z0M, AIRD, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: TAUW, TAUWDIR, PHIWA
      LOGICAL, INTENT(IN) :: LLPHIWA


      INTEGER(KIND=JWIM) :: IJ, M, K, I, J, II

      REAL(KIND=JWRB) :: TAUTOUS2
      REAL(KIND=JWRB) :: COSW, FCOSW2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUHF, PHIHF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: US2, TAUX, TAUY, TAUPX, TAUPY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: USDIRP, UST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SUMT, SUMX, SUMY

      LOGICAL :: LTAUWSHELTER

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)

      DO IJ=KIJS,KIJL
        PHIWA(IJ)   = 0.0_JWRB
        XSTRESS(IJ) = 0.0_JWRB
        YSTRESS(IJ) = 0.0_JWRB
      ENDDO

!*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
!     ------------------------------------------------------------------------

      IF ( LLPHIWA ) THEN
!     full energy flux due to negative Sinput (SL-SPOS)
!     we assume that above NFRE, the contibutions can be neglected
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              PHIWA(IJ) = PHIWA(IJ) + (SL(IJ,K,M)-SPOS(IJ,K,M))*RHOWG_DFIM(M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
!     ---------------------------------------------------------------------------------
      DO M=1,MAXVAL(MIJ(:))
!     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
        K=1
        DO IJ=KIJS,KIJL
          SUMX(IJ) = SPOS(IJ,K,M)*SINTH(K)
          SUMY(IJ) = SPOS(IJ,K,M)*COSTH(K)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            SUMX(IJ) = SUMX(IJ) + SPOS(IJ,K,M)*SINTH(K)
            SUMY(IJ) = SUMY(IJ) + SPOS(IJ,K,M)*COSTH(K)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          CMRHOWGDFTH(IJ) = RHOWGDFTH(IJ,M)*CINV(IJ,M)
          XSTRESS(IJ) = XSTRESS(IJ) + CMRHOWGDFTH(IJ)*SUMX(IJ)
          YSTRESS(IJ) = YSTRESS(IJ) + CMRHOWGDFTH(IJ)*SUMY(IJ)
        ENDDO
      ENDDO

!     TAUW is the kinematic wave stress !
      DO IJ=KIJS,KIJL
        XSTRESS(IJ) = XSTRESS(IJ)/MAX(AIRD(IJ), 1.0_JWRB)
        YSTRESS(IJ) = YSTRESS(IJ)/MAX(AIRD(IJ), 1.0_JWRB)
      ENDDO

      IF ( LLPHIWA ) THEN
        DO M=1,MAXVAL(MIJ(:))
!       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
          K=1
          DO IJ=KIJS,KIJL
            SUMT(IJ) = SPOS(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=KIJS,KIJL
              SUMT(IJ) = SUMT(IJ) + SPOS(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            PHIWA(IJ) = PHIWA(IJ) + RHOWGDFTH(IJ,M)*SUMT(IJ)
          ENDDO
        ENDDO
      ENDIF

!*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
!     ----------------------------------------------------------------------------------

      DO IJ=KIJS,KIJL
        US2(IJ)=UFRIC(IJ)**2
      ENDDO

      IF ( IPHYS == 0 .OR. TAUWSHELTER == 0.0_JWRB) THEN
        LTAUWSHELTER = .FALSE.
        DO IJ=KIJS,KIJL
          USDIRP(IJ)=WDWAVE(IJ)
          UST(IJ)=UFRIC(IJ)
        ENDDO
      ELSE
        LTAUWSHELTER = .TRUE.
        DO IJ=KIJS,KIJL
          TAUX(IJ)=US2(IJ)*SIN(WDWAVE(IJ))
          TAUY(IJ)=US2(IJ)*COS(WDWAVE(IJ))
          TAUPX(IJ)=TAUX(IJ)-TAUWSHELTER*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-TAUWSHELTER*YSTRESS(IJ)
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
          UST(IJ)=(TAUPX(IJ)**2+TAUPY(IJ)**2)**0.25_JWRB
        ENDDO
      ENDIF

      CALL TAU_PHI_HF(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, &
     &                FL1, AIRD, RNFAC,                          &
     &                COSWDIF, SINWDIF2,                         &
     &                UST, TAUHF, PHIHF, LLPHIWA)

      DO IJ=KIJS,KIJL
        XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
        TAUW(IJ) = MAX(TAUW(IJ),0.0_JWRB)
        TAUWDIR(IJ) = ATAN2(XSTRESS(IJ),YSTRESS(IJ))
      ENDDO

      IF ( .NOT. LLGCBZ0) THEN
        TAUTOUS2 = 1.0_JWRB/(1.0_JWRB+EPS1)
        DO IJ=KIJS,KIJL
          TAUW(IJ) = MIN(TAUW(IJ),US2(IJ)*TAUTOUS2)
        ENDDO
      ENDIF

      IF ( LLPHIWA ) THEN
        DO IJ=KIJS,KIJL
          PHIWA(IJ) = PHIWA(IJ) + PHIHF(IJ)
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      END SUBROUTINE STRESSO
