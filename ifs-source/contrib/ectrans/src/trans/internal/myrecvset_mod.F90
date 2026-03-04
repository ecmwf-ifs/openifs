! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MYRECVSET_MOD
CONTAINS
FUNCTION MYRECVSET(KSETS,KMYSET,KSET)


!**** *MYRECVSET* RETURNS SET NUMBER TO SEND TO

!     Purpose.
!     --------
!

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

!

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

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE
INTEGER(KIND=JPIM) :: MYRECVSET
INTEGER(KIND=JPIM),INTENT(IN)  :: KSETS,KMYSET,KSET


!     ------------------------------------------------------------------

!*       1.    Check input argument for validity
!              ---------------------------------

IF(KSETS < 1 .OR. KMYSET > KSETS .OR. KSET > KSETS-1) THEN

  CALL ABORT_TRANS(' MYRECVSET: INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  MYRECVSET = MOD(-KSET-1+KMYSET+KSETS,KSETS)+1

ENDIF

END FUNCTION MYRECVSET
END MODULE MYRECVSET_MOD
