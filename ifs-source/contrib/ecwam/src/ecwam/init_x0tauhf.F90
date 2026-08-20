! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INIT_X0TAUHF

! ----------------------------------------------------------------------

!**** *INIT_X0TAUHF* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR TAU_PHI_HF


!**   INTERFACE.
!     ----------

!       *CALL* *INIT_X0TAUHF

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : X0TAUHF, JTOT_TAUHF, WTAUHF, LLCAPCHNK, LLGCBZ0, LLNORMAGAM
      USE YOWFRED  , ONLY : DELTH
      USE YOWPHYS  , ONLY : BETAMAX  ,ZALP     ,ALPHA, ALPHAMIN, XKAPPA, &
     &                      BETAMAXOXKAPPA2,                             &
     &                      BMAXOKAP, BMAXOKAPDTH, GAMNCONST,            &
     &                      DELTA_THETA_RN
      USE YOWPCONS , ONLY : GM1      ,ZPI

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: J

      REAL(KIND=JWRB) :: CONST1, X0, FF, F, DF, ALPH
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INIT_X0TAUHF',0,ZHOOK_HANDLE)

      BETAMAXOXKAPPA2 = BETAMAX / XKAPPA**2

      BMAXOKAP = DELTA_THETA_RN * BETAMAXOXKAPPA2 /XKAPPA
      BMAXOKAPDTH = BMAXOKAP * DELTH
      GAMNCONST = BMAXOKAP * 0.5_JWRB * ZPI**4 * GM1**3


!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------
      IF(LLGCBZ0 .OR. LLCAPCHNK .OR. LLNORMAGAM ) THEN
        ALPH = ALPHAMIN
      ELSE
        ALPH = ALPHA
      ENDIF

!     find lowest limit for integration X0 *(G/USTAR)
!     ALPHA*X0**2*EXP(XKAPPA/(X0+ZALP))=1
      X0=0.005_JWRB
      DO J=1,30
        FF=EXP(XKAPPA/(X0+ZALP))
        F=ALPH*X0**2*FF-1.0_JWRB
        IF (F.EQ.0.0_JWRB) EXIT
        DF=ALPH*FF*(2.0_JWRB*X0-XKAPPA*(X0/(X0+ZALP))**2)
        X0=X0-F/DF
      ENDDO
      X0TAUHF=X0

      CONST1 = BETAMAXOXKAPPA2/3.0_JWRB

      ! Simpson Integration weights (JTOT_TAUHF must be odd) !
      WTAUHF(1)=CONST1
      DO J=2,JTOT_TAUHF-1,2
        WTAUHF(J)=4.0_JWRB*CONST1
        WTAUHF(J+1)=2.0_JWRB*CONST1
      ENDDO
      WTAUHF(JTOT_TAUHF)=CONST1

      IF (LHOOK) CALL DR_HOOK('INIT_X0TAUHF',1,ZHOOK_HANDLE)

      END SUBROUTINE INIT_X0TAUHF
