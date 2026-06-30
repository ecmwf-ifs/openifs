SUBROUTINE CNT21S
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *CNT21S*  -  Controls integration job at level 2.

!     Purpose.
!     --------
!           Controls integration job at level 2.

!**   Interface.
!     ----------
!        *CALL* *CNT21S

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.   CNT31S - controls integration on level 3
!     ---------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the 1D surface model

!     Author.
!     -------
!        Jean-Francois Mahfouf  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-03-01

!     ------------------------------------------------------------------

IMPLICIT NONE

#include "cnt31s.intfb.h"

!*       1.    CALL LEVEL 3 CONTROL ROUTINE.
!              -----------------------------

CALL CNT31S

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE CNT21S


