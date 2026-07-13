SUBROUTINE CNT31S
USE PARKIND1  ,ONLY : JPIM     ,JPRB,  JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN1S , ONLY : NULOUT
USE YOMDPHY  , ONLY : NPDONE

IMPLICIT NONE
#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *CNT31S*  - Controls integration job at level 3

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *CNT31S

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      Calls SUINIF1S, CNT4
!      Called by CNT2

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the one-column surface IFS

!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-06-16
!     ------------------------------------------------------------
#endif

LOGICAL LOPLEFT   ! More gridpoints to do?

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "suinif1s.intfb.h"
#include "cnt41s.intfb.h"

IF (LHOOK) CALL DR_HOOK('CNT31S',0,ZHOOK_HANDLE)

LOPLEFT=.TRUE.
NPDONE=0
!     ------------------------------------------------------------

!*       1.    Initialize LEVEL 3 COMMONS.
!              ---------------------------

DO 
  WRITE(UNIT=NULOUT,FMT='(A)')' START CNT31S'
  CALL SUINIF1S(LOPLEFT)
!      -----------------------------------------------------------
 
!*       2.    INTEGRATION.
!              ------------

  CALL CNT41S
!      -----------------------------------------------------------

  WRITE(UNIT=NULOUT,FMT='('' END CNT31S'')')
  IF(.NOT.LOPLEFT)EXIT
ENDDO
!     ------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CNT31S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE CNT31S
