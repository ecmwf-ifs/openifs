! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MYSENDSET_MOD
CONTAINS
FUNCTION MYSENDSET(KSETS,KMYSET,KSET)


!**** *MYSENDSET* RETURNS SET NUMBER TO SEND TO

!     Purpose.
!     --------
!

!**   Interface.
!     ----------
!        ISENDSET = MYSENDSET(KSETS,KMYSET,KSET)

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

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE
INTEGER(KIND=JPIM) :: MYSENDSET
INTEGER(KIND=JPIM),INTENT(IN)  :: KSETS,KMYSET,KSET


!     ------------------------------------------------------------------

!*       1.    Check input argument for validity
!              ---------------------------------

IF(KSETS < 1 .OR. KMYSET > KSETS .OR. KSET > KSETS-1) THEN

  CALL ABORT_TRANS(' MYSENDSET: INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  MYSENDSET = MOD(KMYSET+KSET-1,KSETS)+1

ENDIF

END FUNCTION MYSENDSET
END MODULE MYSENDSET_MOD
