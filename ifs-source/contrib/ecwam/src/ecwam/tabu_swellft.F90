! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE TABU_SWELLFT

!**** *TABU_SWELLFT* - FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!     FABRICE ARDHUIN  IFREMER  2013

!*    PURPOSE.
!     --------
!     TO ESTIMATE FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!**   INTERFACE.
!     ----------

!       *CALL*TABU_SWELLFT*

!     METHOD.
!     -------
!       TABULATION ON KELVIN FUNCTIONS.

!     EXTERNALS.
!     -----------

!     KERKEI  (zeroth order Kelvin function Ker and Kei)

! ----------------------------------------------------------------------
!
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWTABL  , ONLY : IAB      ,SWELLFT

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "kerkei.intfb.h"

      INTEGER(KIND=JWIM), PARAMETER :: NITER=100
      REAL(KIND=JWRB),    PARAMETER :: ABMIN=0.3_JWRB
      REAL(KIND=JWRB),    PARAMETER :: ABMAX=8.0_JWRB, KAPPA=0.40_JWRB
!     VARIABLE.   TYPE.     PURPOSE.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
! ----------------------------------------------------------------------
      INTEGER(KIND=JWIM) :: I,ITER
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELAB
      REAL(KIND=JWRB) :: KER, KEI
      REAL(KIND=JWRB) :: ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,DZETA0,DZETA0MEMO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TABU_SWELLFT',0,ZHOOK_HANDLE)

!AR:DZETA0 was not set
      DZETA0 = 0.0_JWRB
!
      DELAB   = (ABMAX-ABMIN)/REAL(IAB,JWRB)
      L10=LOG(10.0_JWRB)
      DO I=1,IAB
         ABRLOG=ABMIN+REAL(I,JWRB)*DELAB
         ABR=EXP(ABRLOG*L10)
         FACT=1/ABR/(21.2_JWRB*KAPPA)
         FSUBW=0.05_JWRB
         DO ITER=1,NITER
            FSUBWMEMO=FSUBW
            DZETA0MEMO=DZETA0
            DZETA0=FACT*FSUBW**(-0.5_JWRB)
            CALL KERKEI(2.0_JWRB*SQRT(DZETA0),KER,KEI)
            FSUBW=0.08_JWRB/(KER**2+KEI**2)
            FSUBW=0.5_JWRB*(FSUBWMEMO+FSUBW)
            DZETA0=0.5_JWRB*(DZETA0MEMO+DZETA0)
         ENDDO   
         SWELLFT(I)=FSUBW
      ENDDO

      IF (LHOOK) CALL DR_HOOK('TABU_SWELLFT',1,ZHOOK_HANDLE)

      END SUBROUTINE TABU_SWELLFT
