! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE TRANS_END(CDMODE)

!**** *TRANS_END* - Terminate transform package 

!     Purpose.
!     --------
!     Terminate transform package. Release all allocated arrays.

!**   Interface.
!     ----------
!     CALL TRANS_END

!     Explicit arguments : None
!     -------------------- 

!     Method.
!     -------

!     Externals.  None
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03


!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
CHARACTER*5, OPTIONAL, INTENT(IN) :: CDMODE


END SUBROUTINE TRANS_END
END INTERFACE
