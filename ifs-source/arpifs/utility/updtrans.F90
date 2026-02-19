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

SUBROUTINE UPDTRANS(KRESOL,LDETRANS)

!**** *UPDTRANS*  - Update transforms setup

!     Purpose.
!     --------
!        To get the LAM/global aspect of a given resolution

!**   Interface.
!     ----------
!        *CALL* *UPDTRANS*

!        Explicit arguments :
!        --------------------

!            KRESOL  : resolution tag (input)
!            LDETRANS: .TRUE. is the geometry corresponding to KRESOL is LAM

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        The spectral transforms package swaps the resolution when invoked with
!        an explicit resolution tag, which causes de-synchronization of pointers
!        when going from Global to LAM, or from LAM to LAM if trans_inq only is
!        used.
!        Calling etrans_inq  systematically would solve the problem but this is
!        not allowed in Ifs ; therefore this subroutine makes the reverse
!        pointers swapping

!     Externals.
!     ----------
!           Spectral transforms subroutines

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 24-Aug-2012

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM    ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KRESOL
LOGICAL           ,INTENT(OUT)   :: LDETRANS

INTEGER(KIND=JPIM) :: ICUR_RESOL
LOGICAL :: LLCUR_LAM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "get_current.h"
#include "trans_inq.h"
#include "etrans_inq.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDTRANS',0,ZHOOK_HANDLE)

CALL GET_CURRENT(KRESOL=ICUR_RESOL,LDLAM=LLCUR_LAM)
IF (LLCUR_LAM) THEN
! Get the target geometry
  CALL ETRANS_INQ(KRESOL=KRESOL,LDLAM=LDETRANS)
! Reverse the resolution swapping:
  CALL ETRANS_INQ(KRESOL=ICUR_RESOL)
ELSE
! Get the target geometry
  CALL TRANS_INQ(KRESOL=KRESOL,LDLAM=LDETRANS)
! Reverse the resolution swapping:
  CALL TRANS_INQ(KRESOL=ICUR_RESOL)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDTRANS',1,ZHOOK_HANDLE)
END SUBROUTINE UPDTRANS
