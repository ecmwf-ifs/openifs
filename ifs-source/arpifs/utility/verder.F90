! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE VERDER(KPROMA,KSTART,KPROF,KFLEV_IN,KFLEV_OUT,PDERI,PIN,POUT)

!**** *VERDER*   Vertical derivative for VFE.

!     Purpose.
!     --------
!          This subroutine computes the vertical derivative (with respect
!          to eta) of a function given at full model
!          levels using a general scheme developed by means of finite-elements

!**   Interface.
!     ----------
!        *CALL* *VERDER(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA   - horizontal dimension.
!          KSTART   - first element of work.
!          KPROF    - depth of work.
!          KFLEV_IN - number of input model layers
!          KFLEV_OUT- number of output model layers
!          PDERI    - derivative operator (first or second order)
!          PIN      - Input field

!        OUTPUT:
!          POUT     - d(PIN)/d(eta) or d2(PIN)/d(eta)2  at each model level

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none

!     Reference.

!     Author.
!     -------
!         M. Hortal (ECMWF) (after VERINT)

!     Modifications.
!     --------------
!         P. Smolikova and J. Vivoda (Oct 2013): more flexible version
!      F. Vana  05-Mar-2015  Support for single precision
!      F. Vana  14-jan-2020  Always use double precision regardless the model precision
!      M. Diamantakis Feb 2021: SP fix for NH runs
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK


!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV_IN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV_OUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
REAL(KIND=JPRD)   ,INTENT(IN)    :: PDERI(KFLEV_OUT,KFLEV_IN)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(KPROMA,KFLEV_IN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(KPROMA,KFLEV_OUT)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDB
INTEGER(KIND=JPIM) ::  JLEV, JROF
REAL(KIND=JPRD) :: ZIN(KPROMA,KFLEV_IN)
REAL(KIND=JPRD) :: ZOUT(KPROMA,KFLEV_OUT)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VERDER',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------


! MATRIX MULTIPLY REPLACED BY BLAS ROUTINE

!DO JLEV=1,KFLEV_OUT
!  DO JROF=KSTART,KPROF
!    POUT(JROF,JLEV)=0._JPRB
!  ENDDO
!ENDDO
!DO JLE=1,KFLEV_IN
!  DO JLEV=1,KFLEV_OUT
!    DO JROF=KSTART,KPROF
!      POUT(JROF,JLEV)=POUT(JROF,JLEV)+PDERI(JLEV,JLE)*PIN(JROF,JLE)
!    ENDDO
!  ENDDO
!ENDDO

IDB=SIZE(PDERI,DIM=1)
IF (JPRB == JPRD) THEN
  ! No conversion required
  CALL DGEMM('N','T',KPROF-KSTART+1,KFLEV_OUT,KFLEV_IN, &
       & 1.0_JPRD,PIN,KPROMA,PDERI,IDB,0.0_JPRD,POUT,KPROMA)  
ELSE
  ! Converting to double
  DO JLEV=1,KFLEV_IN
    DO JROF=KSTART,KPROF
      ZIN(JROF,JLEV)=REAL(PIN(JROF,JLEV),JPRD)
    ENDDO
  ENDDO
  CALL DGEMM('N','T',KPROF-KSTART+1,KFLEV_OUT,KFLEV_IN, &
       & 1.0_JPRD,ZIN,KPROMA,PDERI,IDB,0.0_JPRD,ZOUT,KPROMA)
  DO JLEV=1,KFLEV_OUT
    DO JROF=KSTART,KPROF
      POUT(JROF,JLEV)=REAL(ZOUT(JROF,JLEV),JPRB)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VERDER',1,ZHOOK_HANDLE)
END SUBROUTINE VERDER
