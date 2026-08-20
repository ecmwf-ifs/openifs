! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPABORT(CDMESSAGE) 

! ----------------------------------------------------------------------

!**** *MPABORT* -

!     J. BIDLOT     ECMWF   2001 


!*    PURPOSE.
!     --------
!     CALL ABORT AFTER CLOSING OPENED UNITS 

!**   INTERFACE.
!     ----------
!     *CALL MPABORT(CDMESSAGE)*

!     METHOD.
!     -------
!     CALL MPCLOSE_UNIT THEN CALL ABORT 

!     EXTERNALS.
!     ----------
!       MPCLOSE_UNIT
!       WAM_ABORT

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU  
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpclose_unit.intfb.h"

      CHARACTER(LEN=*),INTENT(IN) :: CDMESSAGE

      LOGICAL :: LLABORT=.TRUE.

      IF (LWCOU) THEN
        CALL WAM_ABORT(CDMESSAGE)
      ELSE
        CALL MPCLOSE_UNIT
        STOP 1
      ENDIF

      END SUBROUTINE MPABORT
