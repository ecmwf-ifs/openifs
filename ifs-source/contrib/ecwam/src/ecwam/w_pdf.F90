REAL(KIND=JWRB) FUNCTION W_PDF(X,AA,BB,C4,C3,LLNOTAIL)
!
!***  DETERMINE ENVELOPE WAVE HEIGHT PROBABILITY DISTRIBUTION            
!
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     X              REAL         NORMALIZED WAVE HEIGHT
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

      REAL(KIND=JWRB), INTENT(IN) :: X,AA,BB,C4,C3
      LOGICAL, INTENT(IN) :: LLNOTAIL

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: AC4,AC3,ARG_A,E,Z,P_E,ZFAC   

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('W_PDF',0,ZHOOK_HANDLE)

!*    1. DETERMINE PDF.
!     -----------------

      IF (LLNOTAIL) THEN
        AC4 = (2._JWRB*X**4-4._JWRB*X**2+1._JWRB)
        AC3 = (4._JWRB*X**6-18._JWRB*X**4+18._JWRB*X**2-3._JWRB)
        ARG_A = EXP(-2._JWRB*X**2)*(1._JWRB+C4*AC4+C3**2*AC3)
        W_PDF = 4._JWRB*X*ARG_A
      ELSE
        E = 2._JWRB*X**2
        Z = -AA+SQRT(AA**2+BB*E)
        P_E = EXP(-Z)
        ZFAC = BB/(2._JWRB*(Z+AA))
        W_PDF = 4._JWRB*X*ZFAC*P_E
      ENDIF

      IF (LHOOK) CALL DR_HOOK('W_PDF',0,ZHOOK_HANDLE)

END FUNCTION W_PDF
