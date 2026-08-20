! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

REAL(KIND=JWRB) FUNCTION TRANSF(XK,D)
!
!***  DETERMINE NARROW BAND LIMIT NONLINEAR TRANSFER FUNCTION            
!     BASED ON TECH MEMO 464 BY P. JANSSEN AND M. ONORATO
!     
!
!     AUTHOR:  P.A.E.M. JANSSEN ECMWF JUNE 2005
!     ------
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     XK             REAL         WAVE NUMBER
!     D              REAL         DEPTH
!
!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G     ,DKMAX

!----------------------------------------------------------------------
 
      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: XK,D

      REAL(KIND=JWRB) :: EPS,X,T_0,OM,C_0,V_G,DV_G,XNL_1,XNL_2,XNL

      EPS=0.0001_JWRB
!
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------
!     
      IF ( D < 999.0_JWRB .AND. D > 0.0_JWRB) THEN
        X   = XK*D
        IF ( X > DKMAX) THEN
          TRANSF = 1.0_JWRB 
        ELSE
          T_0 = TANH(X)
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          IF (X < EPS) THEN
            V_G = 0.5_JWRB*C_0
            V_G = C_0
          ELSE
            V_G = 0.5_JWRB*C_0*(1.0_JWRB+2.0_JWRB*X/SINH(2.0_JWRB*X))
          ENDIF
          DV_G = (T_0-X*(1.0_JWRB-T_0**2))**2+4.0_JWRB*X**2*T_0**2*(1.0_JWRB-T_0**2)
       
          XNL_1 = (9.0_JWRB*T_0**4-10.0_JWRB*T_0**2+9.0_JWRB)/(8.0_JWRB*T_0**3)
          XNL_2 = ((2.0_JWRB*V_G-0.5_JWRB*C_0)**2/(G*D-V_G**2)+1.0_JWRB)/X

          XNL = XNL_1-XNL_2
          TRANSF = XNL**2/(DV_G*T_0**8)
        ENDIF
      ELSE
        TRANSF = 1.0_JWRB 
      ENDIF
!
END FUNCTION TRANSF
