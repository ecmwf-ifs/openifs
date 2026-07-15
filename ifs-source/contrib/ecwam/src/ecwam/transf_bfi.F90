! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

REAL(KIND=JWRB) FUNCTION TRANSF_BFI(XK0,D,XNU,SIG_TH)
 
!***  DETERMINE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR
!     THE FINITE DEPTH CASE             
 
 
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XK0            REAL         WAVE NUMBER
!     D              REAL         DEPTH
 
!     XNU            REAL         RELATIVE SPECTRAL WIDTH
!     SIG_TH         REAL         RELATIVE WIDTH in DIRECTION
 
!     BF**2 = (2 S^2)/SIG_OM^2) . (2 V_G/C_0)^2 . T_0/K_0^3 
 
!----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G     ,DKMAX
      USE YOWSHAL  , ONLY : BATHYMAX, XKDMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
                     
!----------------------------------------------------------------------

      IMPLICIT NONE
      REAL(KIND=JWRB), INTENT(IN) :: XK0,D,XNU,SIG_TH

      REAL(KIND=JWRB), PARAMETER :: EPS=0.0001_JWRB
      REAL(KIND=JWRB), PARAMETER :: TRANSF_BFI_MIN = -4._JWRB
      REAL(KIND=JWRB), PARAMETER :: TRANSF_BFI_MAX = 4._JWRB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X,XK,T_0,T_0_SQ,OM,C_0,C_S_SQ,V_G,V_G_SQ,D2OM
      REAL(KIND=JWRB) :: XNL_1,XNL_2,XNL_3,XNL_4
      REAL(KIND=JWRB) :: XNL,ALP,ZFAC,T_NL
!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TRANSF_BFI',0,ZHOOK_HANDLE)

!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------
     
      IF (D < BATHYMAX .AND. D > 0._JWRB) THEN
        X   = XK0*D
        IF ( X .GT. DKMAX) THEN
          TRANSF_BFI = 1._JWRB 
        ELSE
          XK  = MAX(XK0,XKDMIN/D)
          X   = XK*D
          T_0 = TANH(X)
          T_0_SQ = T_0**2
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          C_S_SQ = G*D
          IF (X < EPS) THEN
            V_G = C_0
          ELSE
            V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          ENDIF
          V_G_SQ = V_G**2
          D2OM = (T_0-X*(1._JWRB-T_0_SQ))**2+4._JWRB*X**2*T_0_SQ*(1._JWRB-T_0_SQ)
      
          XNL_1 = (9._JWRB*T_0_SQ**2-10._JWRB*T_0_SQ+9._JWRB)/(8._JWRB*T_0_SQ*T_0)
          XNL_2 = ((2._JWRB*V_G-0.5_JWRB*C_0)**2/(G*D-V_G_SQ)+1._JWRB)/X
          XNL_4 = 1./(4._JWRB*T_0)*(2._JWRB*C_0+V_G*(1._JWRB-T_0_SQ))**2/(C_S_SQ-V_G_SQ)
          ALP = (1._JWRB-V_G_SQ/C_S_SQ)*C_0**2/V_G_SQ
          ZFAC = SIG_TH**2/(SIG_TH**2+ALP*XNU**2)
          XNL_3 = ZFAC*XNL_4

!!          T_NL = XK**3*(XNL_1-XNL_2+XNL_3)
!!          TRANSF_BFI = 4._JWRB*(V_G/C_0)**2*T_NL/XK**3*T_0/D2OM
          T_NL = XNL_1-XNL_2+XNL_3
          TRANSF_BFI = 4._JWRB*(V_G/C_0)**2*T_NL*T_0/D2OM
          TRANSF_BFI = MAX(MIN(TRANSF_BFI_MAX,TRANSF_BFI),TRANSF_BFI_MIN)
        ENDIF
      ELSE
        TRANSF_BFI = 1._JWRB
      ENDIF

      IF (LHOOK) CALL DR_HOOK('TRANSF_BFI',1,ZHOOK_HANDLE)

END FUNCTION TRANSF_BFI
