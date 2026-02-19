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

FUNCTION MYRECVSET(KSETS,KMYSET,KSET)

!**** *MYRECVSET* RETURNS SET NUMBER TO SEND TO

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        ISENDSET = MYRECVSET(KSETS,KMYSET,KSET)

!        Explicit arguments :  
!        --------------------
!                  input:   KSETS

!        Implicit arguments :  NONE
!        --------------------
!     Method.
!     -------

!     Externals.
!     ----------
!         NONE

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
INTEGER(KIND=JPIM) :: MYRECVSET
INTEGER(KIND=JPIM),INTENT(IN)  :: KSETS,KMYSET,KSET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.    Check input argument for validity 
!              ---------------------------------

IF (LHOOK) CALL DR_HOOK('MYRECVSET',0,ZHOOK_HANDLE)
IF(KSETS < 1 .OR. KMYSET > KSETS .OR. KSET > KSETS-1) THEN

  CALL ABOR1(' MYRECVSET: INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  MYRECVSET = MOD(-KSET-1+KMYSET+KSETS,KSETS)+1

ENDIF

IF (LHOOK) CALL DR_HOOK('MYRECVSET',1,ZHOOK_HANDLE)
END FUNCTION MYRECVSET
