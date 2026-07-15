REAL(KIND=JWRB) FUNCTION  W_CDFMAX(E,AA,BB,N,H_C,C4,C3,LLNOTAIL)
!
!***  DETERMINE MAX WAVE HEIGHT DISTRIBUTION            
!
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     E              REAL         NORMALIZED WAVE ENERGY
!     AA             REAL         RESIDORI FIT PARAMETER
!     BB             REAL         BB = 2*(AA+1)
!     N              INTEGER      NUMBER OF DEGREES OF FREEDOM 
!     C4             REAL         KURTOSIS
!     C3             REAL         SKEWNESS
!     LLNOTAIL         LOGICAL      IF TRUE EDGWORTH DISTRIBUTION IS USED
!                                 OTHERWISE THE STRETCHED EXPONENTIAL 
!                                 DISTRIBUTION IS USED
!
!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: N   
      REAL(KIND=JWRB), INTENT(IN) :: E,AA,BB,H_C,C4,C3
      LOGICAL, INTENT(IN) :: LLNOTAIL

      REAL(KIND=JWRB), PARAMETER :: EPS = 0.0000000001_JWRB
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: H,XN,AE,BE,P_E,Z

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('W_CDFMAX',0,ZHOOK_HANDLE)

!*    1. DETERMINE W_CDFMAX.
!     -------------------
      H =SQRT(0.5_JWRB*E)
      XN = REAL(N)*H

      IF (LLNOTAIL) THEN
        AE = 0.5_JWRB*E*(E-2._JWRB)
        BE = 0.5_JWRB*E*(E**2-6._JWRB*E+6._JWRB)
        P_E = EXP(-E)*(1._JWRB+C4*AE+C3**2*BE)

        W_CDFMAX = 1._JWRB-EXP(-XN*P_E)
      ELSE    
        Z = -AA+SQRT(AA**2+BB*E)
        P_E = EXP(-Z)

        W_CDFMAX = 1._JWRB-EXP(-XN*P_E)

      ENDIF

      W_CDFMAX = MAX(EPS,W_CDFMAX)

      IF (LHOOK) CALL DR_HOOK('W_CDFMAX',1,ZHOOK_HANDLE)

END FUNCTION W_CDFMAX
