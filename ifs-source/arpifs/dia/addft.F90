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

SUBROUTINE ADDFT(KLUNOUT,CDFT,CDVAR,CDNAME,K)
!     Purpose.
!     --------
!     Adding a new flux/tendency description
!     

!   Interface.
!   ----------
!      CALL ADDFT(...)

!    Explicit arguments
!    ------------------
! KLUNOUT  : logical number of standard uotput unit
! CDFT     : "F" if flux, "T" if tendency
! CDVAR    : variable that is changed by this term
! CDNAME   : name of flux/tendency
! K        : ordinal number of this flux/tendency

!     Author.
!     -------
!      Tomislav Kovacic

!     Modifications.
!     --------------
!      Original : 2006-03-23
!      Modifications:
!    -------------------------------------------------------------------------

USE PARKIND1, ONLY  : JPIM, JPRB
USE YOMPHFT,  ONLY  : YAPFT, NDDHFT
USE YOMHOOK,  ONLY  : LHOOK, DR_HOOK, JPHOOK
!    -------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)   :: KLUNOUT
CHARACTER(LEN=1),   INTENT(IN)   :: CDFT
CHARACTER(LEN=2),   INTENT(IN)   :: CDVAR
CHARACTER(LEN=10),  INTENT(IN)   :: CDNAME
INTEGER(KIND=JPIM), INTENT(INOUT):: K

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!    -------------------------------------------------------------------------
#include "abor1.intfb.h"
!    -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ADDFT',0,ZHOOK_HANDLE)

K=K+1
IF ( K>NDDHFT ) THEN
  WRITE(KLUNOUT,'(" NDDHFT= ",I5)') NDDHFT
  CALL ABOR1('ADDFT: NDDHFT is too small!')
ENDIF
YAPFT(K)%CFT  = CDFT
YAPFT(K)%CVAR = CDVAR
YAPFT(K)%CNAME= CDNAME

IF (LHOOK) CALL DR_HOOK('ADDFT',1,ZHOOK_HANDLE)

END SUBROUTINE ADDFT
