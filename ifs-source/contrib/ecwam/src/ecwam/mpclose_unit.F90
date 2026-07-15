! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPCLOSE_UNIT

! ----------------------------------------------------------------------

!**** *MPCLOSE_UNIT* -

!     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING
!                            2007   ONLY USE IT TO FLUSH ALL UNITS.

!*    PURPOSE.
!     --------
!     CLOSE ALL OPENED UNIT

!**   INTERFACE.
!     ----------
!     *CALL MPCLOSE_UNIT*

!     METHOD.
!     -------
!     LOOPS ON ALL UNIT NUMBER AND INQUIRES IF UNIT IS OPENED 
!     AND THEN CLOSES IT

!     EXTERNALS.
!     ----------
!       NONE.
!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IU
      LOGICAL :: LOPENED

! ----------------------------------------------------------------------

      DO IU=1,99
        INQUIRE(UNIT=IU,OPENED=LOPENED)
        IF (LOPENED) THEN 
          CALL FLUSH(IU)
!!!!          CLOSE(IU)
        ENDIF
      ENDDO 
      END SUBROUTINE MPCLOSE_UNIT
