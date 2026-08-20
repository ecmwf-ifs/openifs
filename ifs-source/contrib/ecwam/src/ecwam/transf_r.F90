! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

REAL(KIND=JWRB) FUNCTION TRANSF_R(XK0,D)
 
!***  DETERMINE DEPTH-DEPENDENT RATIO OF ANGULAR WIDTH AND FREQUENCY WIDTH
 
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XK0            REAL         WAVE NUMBER
!     D              REAL         DEPTH
 
 
!     FOR R SEE EQ. (A14) OF TECH MEMO 813 
 
!----------------------------------------------------------------------
 
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G     ,DKMAX
      USE YOWSHAL  , ONLY : BATHYMAX, XKDMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), PARAMETER :: EPS=0.0001_JWRB
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: XK0,D
      REAL(KIND=JWRB) :: X,XK,T_0,T_0_SQ,OM,C_0,V_G,D2OM

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TRANSF_R',0,ZHOOK_HANDLE)
 
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------

      IF (D < BATHYMAX .AND. D > 0._JWRB .AND. XK0 > 0._JWRB) THEN
        X = XK0*D
        IF ( X > DKMAX) THEN
          TRANSF_R = 0.5_JWRB 
        ELSE
          XK  = MAX(XK0,XKDMIN/D)
          X   = XK*D
          T_0 = TANH(X)
          T_0_SQ = T_0**2
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          IF (X < EPS) THEN
            V_G = C_0
          ELSE
            V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          ENDIF
          D2OM = (T_0-X*(1._JWRB-T_0_SQ))**2+4._JWRB*X**2*T_0_SQ*(1._JWRB-T_0_SQ)
          TRANSF_R = 4._JWRB*(V_G/C_0)**3*T_0_SQ/D2OM
        ENDIF
      ELSE
        TRANSF_R = 0.5_JWRB
      ENDIF
 
      IF (LHOOK) CALL DR_HOOK('TRANSF_R',1,ZHOOK_HANDLE)

END FUNCTION TRANSF_R
