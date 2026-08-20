! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDIWBK (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)

! ----------------------------------------------------------------------

!**** *SDIWBK* - COMPUTATION OF BOTTOM-INDUCED WAVE BREAKING DISSIPATION


!*    PURPOSE.
!     --------
!       COMPUTE BOTTOM-INDUCED DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDIWBK (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1*     - SPECTRUM.
!          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*      - TOTAL SOURCE FUNCTION ARRAY
!          *DEPTH*   - WATER DEPTH
!          *EMAXDPT* - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
!          *EMEAN*   - MEAN ENERGY DENSITY
!          *F1MEAN*  - MEAN FREQUENCY BASED ON 1st MOMENT.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE,    NFRE_RED
      USE YOWSTAT  , ONLY : LBIWBK

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL),INTENT(IN):: DEPTH, EMAXDPT, EMEAN, F1MEAN


      INTEGER(KIND=JWIM) :: IJ, K, M, IC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ALPH, ARG, Q, Q_OLD, REL_ERR, EXPQ
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: SDS

      REAL, PARAMETER :: ALPH_B_J = 1.0_JWRB
      REAL, PARAMETER :: COEF_B_J=2*ALPH_B_J
      REAL, PARAMETER :: DEPTHTRS=50.0_JWRB 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDIWBK',0,ZHOOK_HANDLE)

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      IF (LBIWBK) THEN
!       (FOLLOWING BATTJES-JANSSEN AND BEJI)
        DO IJ=KIJS,KIJL
           IF(DEPTH(IJ) < DEPTHTRS) THEN
             ALPH = 2.0_JWRB*EMAXDPT(IJ)/EMEAN(IJ)
             ARG  = MIN(ALPH,50.0_JWRB)
             Q_OLD = EXP(-ARG)
!            USE NEWTON-RAPHSON METHOD
             DO IC=1,15
               EXPQ = EXP(-ARG*(1.0_JWRB-Q_OLD))
               Q = Q_OLD - (EXPQ-Q_OLD)/(ARG*EXPQ-1.0_JWRB)
               REL_ERR=ABS(Q-Q_OLD)/Q_OLD
               IF(REL_ERR.LT.0.00001_JWRB) EXIT
               Q_OLD = Q
             ENDDO
             Q=MIN(Q,1.0_JWRB)
             SDS(IJ) = COEF_B_J*ALPH*Q*F1MEAN(IJ)
           ENDIF
        ENDDO 
      
        DO M = 1, NFRE_RED
           DO K=1,NANG
              DO IJ=KIJS,KIJL
                IF(DEPTH(IJ) < DEPTHTRS) THEN
                  SL(IJ,K,M) = SL(IJ,K,M)-SDS(IJ)*FL1(IJ,K,M)
                  FLD(IJ,K,M) = FLD(IJ,K,M)-SDS(IJ)
                ENDIF
              ENDDO
           ENDDO
        ENDDO
     
      ENDIF

      IF (LHOOK) CALL DR_HOOK('SDIWBK',1,ZHOOK_HANDLE)

      END SUBROUTINE SDIWBK
