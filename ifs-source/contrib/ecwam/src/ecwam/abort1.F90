! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE ABORT1

! ----------------------------------------------------------------------

!**** *ABORT1* -

!     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING


!*    PURPOSE.
!     --------
!     TO ABORT AFTER FLUSHING OPENED UNITS

!**   INTERFACE.
!     ----------
!     *CALL ABORT1*

!     METHOD.
!     -------
!     CALL WAM_ABORT

!     EXTERNALS.
!     ----------
!       WAM_ABORT

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWTEST, ONLY: IU06
      USE YOWABORT, ONLY : WAM_ABORT
! ----------------------------------------------------------------------

      IMPLICIT NONE

      CALL FLUSH (IU06)
      CALL WAM_ABORT()

      END SUBROUTINE ABORT1
