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

SUBROUTINE SPNORMBE(YDLAP,YDLEP,YDDIM,YDEDIM,PX,PY,PSN,PSM,KLEV,PMET,KFLEV,KFLSUR,LDNWAVE)

!**** *SPNORMBE* - Compute norms in spectral space for KE, multitasking interface

!     Purpose.
!     --------
!        Compute the norm in spectral space, used for KE

!**   Interface.
!     ----------
!        *CALL* *SPNORMBE(...)

!        Explicit arguments :
!        --------------------

!        PX      : Input array nr 1
!        PY      : Input array nr 2
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
!      Calls SPNORMBM.
!      Called by SPNORM.

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
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(TLAP)        , INTENT(IN)    :: YDLAP
TYPE(TLEP)        , INTENT(IN)    :: YDLEP
TYPE(TDIM)        , INTENT(IN)    :: YDDIM
TYPE(TEDIM)       , INTENT(IN)    :: YDEDIM
INTEGER(KIND=JPIM), INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM), INTENT(IN)    :: KFLSUR 
REAL(KIND=JPRB)   , INTENT(IN)    :: PX(:,:) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PY(:,:) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSN(:,:) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSM(:,:) 
INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   , INTENT(IN)    :: PMET(:) 
LOGICAL           , INTENT(IN)    :: LDNWAVE 
REAL(KIND=JPRB) :: ZSN(KFLEV,YDDIM%NUMP), ZSM(KFLEV,YDDIM%NUMP)

INTEGER(KIND=JPIM) :: JLEV, JML
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "spnormbm.intfb.h"

!     ------------------------------------------------------------------

!*       1.    Call SPNORMBM.
!              --------------

IF (LHOOK) CALL DR_HOOK('SPNORMBE',0,ZHOOK_HANDLE)
ASSOCIATE(NUMP=>YDDIM%NUMP)
CALL SPNORMBM(YDLAP,YDLEP,YDDIM,YDEDIM,PX,PSN,PSM,KLEV,PMET,KFLEV,KFLSUR,LDNWAVE)

CALL SPNORMBM(YDLAP,YDLEP,YDDIM,YDEDIM,PY,ZSN,ZSM,KLEV,PMET,KFLEV,KFLSUR,LDNWAVE)

IF (LDNWAVE) THEN
  DO JLEV=1,KLEV
    DO JML=1,NUMP
      PSN(JLEV,JML)=(PSN(JLEV,JML)+ZSN(JLEV,JML))*0.5_JPRB
    ENDDO
  ENDDO
ENDIF
DO JLEV=1,KLEV
  DO JML=1,NUMP
    PSM(JLEV,JML)=(PSM(JLEV,JML)+ZSM(JLEV,JML))*0.5_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPNORMBE',1,ZHOOK_HANDLE)
END SUBROUTINE SPNORMBE
