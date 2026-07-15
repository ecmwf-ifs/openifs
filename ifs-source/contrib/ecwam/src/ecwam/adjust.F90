! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE ADJUST (WEST, EAST)

! ----------------------------------------------------------------------

!**** *ADJUST* - ROUTINE TO CORRECT BORDERS OF INTERVALS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.

!**   INTERFACE.
!     ----------

!       *CALL* *ADJUST (WEST, EAST)*
!          *WEST*    - LEFT INTERVAL BORDER IN DEGREE.
!          *EAST*    - RIGHT INTERVAL BORDER IN DEGREE.

!     METHOD.
!     -------

!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT (INOUT) :: WEST, EAST
      REAL(KIND=JWRB) :: OLD_WEST, OLD_EAST

!* 1. CORRECT BORDERS.
!     ----------------

      OLD_WEST = WEST
      OLD_EAST = EAST
      WEST = MOD (WEST + 720._JWRB, 360._JWRB)
      EAST = MOD (EAST + 720._JWRB, 360._JWRB)
      IF (OLD_WEST.NE.OLD_EAST) THEN
        IF (WEST.GE.EAST) WEST = WEST - 360._JWRB
      ENDIF

      END SUBROUTINE ADJUST
