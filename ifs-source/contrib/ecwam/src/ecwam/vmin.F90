! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

REAL(KIND=JWRB) FUNCTION VMIN(XI,XJ,XK,THI,THJ,THK)
!
!***  *VMIN*  DETERMINES THE SECOND-ORDER TRANSFER COEFFICIENT FOR 
!             THREE WAVE INTERACTIONS OF GRAVITY WAVES.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)
!
!     INTERFACE.
!     ----------
!              *VMIN(XI,XJ,XK)*
!                      *XI*   - WAVE NUMBER
!                      *XJ*   - WAVE NUMBER
!                      *XK*   - WAVE NUMBER
!                      *THI*  - WAVE DIRECTION
!                      *THJ*  - WAVE DIRECTION
!                      *THK*  - WAVE DIRECTION
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS, ONLY : G

!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: XI, XJ, XK, THI, THJ, THK

      REAL(KIND=JWRB) :: DEL1, RI, RJ, RK,                              &
     &                   OI, OJ, OK, QI, QJ, QK, RIJ,                   &
     &                   RIK, RJK, SQIJK, SQIKJ, SQJKI, ZCONST, OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL1 = 10.0_JWRB**(-12)
      ZCONST=1.0_JWRB/(4.0_JWRB*SQRT(2.0_JWRB))

      RI = XI
      RJ = XJ
      RK = XK

      OI=OMEG(RI)+DEL1
      OJ=OMEG(RJ)+DEL1
      OK=OMEG(RK)+DEL1

      QI=OI**2/G
      QJ=OJ**2/G
      QK=OK**2/G

      RIJ = RI*RJ*COS(THJ-THI)
      RIK = RI*RK*COS(THK-THI)
      RJK = RJ*RK*COS(THK-THJ)

      SQIJK=SQRT(G*OK/(OI*OJ))
      SQIKJ=SQRT(G*OJ/(OI*OK))
      SQJKI=SQRT(G*OI/(OJ*OK))

      VMIN=ZCONST*( (RIJ-QI*QJ)*SQIJK + (RIK-QI*QK)*SQIKJ + (RJK+QJ*QK)*SQJKI )

END FUNCTION VMIN
