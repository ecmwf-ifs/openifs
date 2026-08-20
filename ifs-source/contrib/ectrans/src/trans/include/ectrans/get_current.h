! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

INTERFACE
SUBROUTINE GET_CURRENT(KRESOL,LDLAM)

!**** *GET_CURRENT* - Extract current information from the transform package

!     Purpose.
!     --------
!     Interface routine for extracting current information from the T.P.

!**   Interface.
!     ----------
!     CALL GET_CURRENT(...)

!     Explicit arguments : (all optional)
!     --------------------
!     KRESOL   - Current resolution
!     LDLAM    -  .T. if the corresponding resolution is LAM, .F. if it is global

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        Ryad El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original : 24-Aug-2012

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM)  ,OPTIONAL,INTENT(OUT)  :: KRESOL
LOGICAL             ,OPTIONAL,INTENT(OUT)  :: LDLAM

END SUBROUTINE GET_CURRENT
END INTERFACE
