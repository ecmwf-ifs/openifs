REAL(KIND=JWRB) FUNCTION W_PMAX(X,AA,BB,N,C4,C3,LLNOTAIL)

!***  DETERMINE MAX WAVE HEIGHT DISTRIBUTION            


!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------

!     X              REAL         NORMALIZED ENVELOPE WAVE HEIGHT
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
      REAL(KIND=JWRB), INTENT(IN) :: X,AA,BB,C4,C3
      LOGICAL, INTENT(IN) :: LLNOTAIL

      REAL(KIND=JWRB), PARAMETER :: EPS = 0.0000000001_JWRB
      REAL(KIND=JWRB) :: XN,AC4,AC3,BC4,BC3,ARG_A,ARG_B
      REAL(KIND=JWRB) :: E,Z,P_E,Y,ZFAC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('W_PMAX',0,ZHOOK_HANDLE)

!*    1. DETERMINE W_PMAX.
!     -----------------
  
      XN = REAL(N)

      IF (LLNOTAIL) THEN
        AC4 = (2._JWRB*X**4-4._JWRB*X**2+1._JWRB)
        AC3 = (4._JWRB*X**6-18._JWRB*X**4+18._JWRB*X**2-3._JWRB)
        BC4 = 2._JWRB*X**2*(X**2-1._JWRB)
        BC3 = 2._JWRB*X**2*(2._JWRB*X**4-6._JWRB*X**2+3._JWRB)
        ARG_A = EXP(-2._JWRB*X**2)*(1._JWRB+C4*AC4+C3**2*AC3)
        ARG_B = EXP(-2._JWRB*X**2)*(1._JWRB+C4*BC4+C3**2*BC3)

        W_PMAX = 4._JWRB*XN*X*ARG_A*EXP(-XN*ARG_B)
      ELSE    
        E = 2._JWRB*X**2
        Z = -AA+SQRT(AA**2+BB*E)
        P_E = EXP(-Z)
        Y = XN*X*P_E
        ZFAC = MAX(0._JWRB,BB/(Z+AA)-1._JWRB/E)
        W_PMAX = 2._JWRB*X*ZFAC*Y*EXP(-Y)
      ENDIF

      W_PMAX = MAX(EPS,W_PMAX)

      IF (LHOOK) CALL DR_HOOK('W_PMAX',1,ZHOOK_HANDLE)

END FUNCTION W_PMAX
