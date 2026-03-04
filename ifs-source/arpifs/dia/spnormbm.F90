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

SUBROUTINE SPNORMBM(YDLAP,YDLEP,YDDIM,YDEDIM,PX,PSN,PSM,KLEV,PMET,KFLEV,KFLSUR,LDNWAVE)

!**** *SPNORMBM* - Compute norms in spectral space, arpifs/aladin interface

!     Purpose.
!     --------
!        Compute the norm in spectral space

!**   Interface.
!     ----------
!        *CALL* *SPNORMBM(...)

!        Explicit arguments :
!        --------------------
!        PX      : Input array
!        PSN     : spectrum at constant n
!        PSM     : spectrum at constant m
!        KLEV    : number of levels of computation
!        PMET    : Metric
!        KFLEV   : first dimensioning of output arrays.
!        KFLSUR  : first dimensioning of input array.
!        LDNWAVE : .TRUE. to compute spectrum at constant n

!        Implicit arguments :  none.
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!      Calls SPNORMB.
!      Called by SPNORMBE and SPNORMBL.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Philippe Courtier  *ECMWF*
!      Original : 92-12-20

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE YOMLAP   , ONLY : TLAP
USE YEMLAP   , ONLY : TLEP
USE YOMDIM   , ONLY : TDIM
USE YEMDIM   , ONLY : TEDIM
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM

IMPLICIT NONE

TYPE(TLAP)        , INTENT(IN)  :: YDLAP
TYPE(TLEP)        , INTENT(IN)  :: YDLEP
TYPE(TDIM)        , INTENT(IN)  :: YDDIM
TYPE(TEDIM)       , INTENT(IN)  :: YDEDIM
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV 
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLSUR 
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEV 
REAL(KIND=JPRB)   , INTENT(IN)  :: PX(:,:) 
REAL(KIND=JPRB)   , INTENT(OUT) :: PSN(:,:) 
REAL(KIND=JPRB)   , INTENT(OUT) :: PSM(:,:) 
REAL(KIND=JPRB)   , INTENT(IN)  :: PMET(:) 
LOGICAL           , INTENT(IN)  :: LDNWAVE 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "espnormb.intfb.h"
#include "spnormb.intfb.h"

!     ------------------------------------------------------------------

!*       1.    CALL SPNORMB.
!              -------------

IF (LHOOK) CALL DR_HOOK('SPNORMBM',0,ZHOOK_HANDLE)

CALL GSTATS(1048,0)

IF (.NOT.LELAM) THEN
  CALL SPNORMB(YDLAP,YDDIM,PX,KLEV,PSN,PSM,PMET,KFLEV,KFLSUR,LDNWAVE)
ELSE
  CALL ESPNORMB(YDLAP,YDLEP,YDDIM,YDEDIM,PX,KLEV,PSN,PSM,PMET,KFLEV,KFLSUR,LDNWAVE)
ENDIF

CALL GSTATS(1048,1)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPNORMBM',1,ZHOOK_HANDLE)
END SUBROUTINE SPNORMBM
