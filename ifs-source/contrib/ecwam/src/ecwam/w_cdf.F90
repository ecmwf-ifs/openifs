REAL(KIND=JWRB) FUNCTION W_CDF(E,AA,BB,C4,C3,LLNOTAIL)
!
!***  DETERMINE ENVELOPE WAVE HEIGHT PROBABILITY            
!
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     E              REAL         NORMALIZED WAVE ENERGY
!     AA             REAL         RESIDORI FIT PARAMETER
!     BB             REAL         BB = 2*(AA+1)
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

      REAL(KIND=JWRB), INTENT(IN) :: E,AA,BB,C4,C3
      LOGICAL, INTENT(IN) :: LLNOTAIL

      INTEGER(KIND=JWIM) :: N   
      REAL(KIND=JWRB), PARAMETER :: EPS = 0.0000000001_JWRB
      REAL(KIND=JWRB) :: H,AE,BE,P_E,Z
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
 
!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('W_CDF',0,ZHOOK_HANDLE)

!*    1. DETERMINE CDF.
!     -----------------

      IF (LLNOTAIL) THEN
        AE = 0.5_JWRB*E*(E-2._JWRB)
        BE = 0.5_JWRB*E*(E**2-6._JWRB*E+6._JWRB)
        W_CDF = EXP(-E)*(1._JWRB+C4*AE+C3**2*BE)
      ELSE    
        Z = -AA+SQRT(AA**2+BB*E)
        W_CDF = EXP(-Z)
      ENDIF

      W_CDF = MAX(EPS,W_CDF)

      IF (LHOOK) CALL DR_HOOK('W_CDF',1,ZHOOK_HANDLE)

END FUNCTION W_CDF
