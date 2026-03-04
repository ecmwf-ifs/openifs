! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!=======================================================================

      SUBROUTINE CONFILE (kuso, kunit, cdate, yafileid, kfail, losuvi)

!-----------------------------------------------------------------------

!****  *confile*  *data aquisition*

!     b. hansen           ECMWF       july    1991

!     PURPOSE.
!     --------
!     Opens data file unformatted for the use of a program.
!     handling of files following the naming convention
!     XXXYYYYMMDDHHMM
!     where XXX is the yafileid
!           YYYY  is the year
!           MM  is the month
!           DD  is the day
!           HH  is the hour
!           MM  is the minute

!**   INTERFACE.
!     ----------
!        CALL  confile (kuso, kunit, cdate, yafileid, kfail)

!        *kuso*       INTEGER       FORTRAN UNIT FOR THE PRINTER OUTPUT.
!        *kunit*      INTEGER       FORTRAN UNIT FOR THE REQUESTED FILE.
!        *cdate*      CHAR*14       DATE TIME GROUP (YYYYMMDDHHMMSS).
!                               !!! The seconds information will not be
!                               !!! used for the filename definition.
!        *yafileid*   CHARACTER*3   FILE ID.
!        *kfail*      INTEGER       = 0 NO ERROR.
!                                   > 0 ERROR.
!        *losuvi*     LOGICAL       IF TRUE THAN MESSAGES SENT TO KUSO

!     EXTERNALS.
!     ----------
!        NONE.

!     METHOD.
!     -------
!        "confile" checks if the required file resides in the working
!         catalog.
!        If not "confile" will return with kfail > 0.
!        If the requested file is not already opened this will be done
!        on unit kunit. IF the unit is already in use it will be 
!        released.

!     REFERENCES.
!     -----------
!        NONE

!----------------------------------------------------------------------

!* 0.    DEFINITIONS.
!  ------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN)  :: kuso, kunit
      INTEGER(KIND=JWIM), INTENT(OUT) :: kfail

      CHARACTER(LEN= 3) :: yafileid
      CHARACTER(LEN=14) :: cdate

      LOGICAL, INTENT(IN) :: losuvi


      INTEGER(KIND=JWIM) :: knum
      CHARACTER(LEN=  7) :: csubna
      CHARACTER(LEN= 12) :: cdate_short
      CHARACTER(LEN= 15) :: yoname
      CHARACTER(LEN=130) :: younitnam

      LOGICAL :: ex, od

      DATA csubna /"confile"/

!----------------------------------------------------------------------

!* 1.0   INITIALIZATION.
!  ---------------------

      kfail = 0

!----------------------------------------------------------------------

!* 2.    CONSTRUCT FILE NAME.
!  --------------------------

      cdate_short=cdate(1:12)
      IF (losuvi) THEN
        WRITE(kuso,*) ' FILE INQUIRY FOR ',yafileid,cdate_short
      ENDIF
      yoname = yafileid
      yoname(4:15) = cdate_short

!----------------------------------------------------------------------

!* 3.0   CHECK STATUS OF UNIT.
!  ---------------------------

      INQUIRE ( UNIT=kunit, OPENED=od, NAME=younitnam, ERR=8200)
      IF (od) THEN
        IF (len_trim(younitnam) .EQ. 15 ) THEN
          IF ( younitnam(1:15) .EQ. yoname ) THEN
            REWIND kunit
            RETURN
          ENDIF
        ELSE
          CLOSE (UNIT=kunit, ERR=8300)
          IF (losuvi) THEN
            WRITE(kuso,*) ' UNIT ',kunit,' WAS OPENED '
            WRITE(kuso,*) ' AND CONNECTED WITH FILE ',younitnam(1:50)
            WRITE(kuso,*) ' IS NOW CLOSED '
          ENDIF
        ENDIF
      ENDIF

!----------------------------------------------------------------------

!* 4.0   CHECK STATUS OF FILE.
!  ---------------------------

      INQUIRE (FILE=yoname, EXIST=ex, OPENED=od, NUMBER=knum, ERR =8100)
      IF (.NOT.ex) THEN
        kfail = 1
        RETURN
      ELSE
        OPEN (FILE=yoname, UNIT=kunit, FORM='unformatted', ERR=8400)
        IF (losuvi) THEN
          WRITE(kuso,*) ' FILE ',yoname,                                &
     &                  ' IS NOW OPENED ON UNIT ', kunit
          WRITE(kuso,*) '  '
        ENDIF
      ENDIF
      RETURN

!----------------------------------------------------------------------

!* 8.0   RETURN IN CASE OF FAILURE.
!  --------------------------------
 8100 CONTINUE
      kfail=20
      WRITE(kuso,*) ' INQUIRY ON FILE NAMED ',yoname,' FAILED'
      WRITE(kuso,*) '  '
      RETURN
 8200 CONTINUE
      kfail=30
      WRITE(kuso,*) ' INQUIRY ON UNIT ',kunit,' FAILED'
      WRITE(kuso,*) '  '
      RETURN
 8300 CONTINUE
      kfail=40
      WRITE(kuso,*) ' CLOSE ON UNIT ',kunit,' FAILED'
      WRITE(kuso,*) '  '
      RETURN
 8400 CONTINUE
      kfail=50
      WRITE(kuso,*) ' OPEN ON FILE NAMED ',yoname,' FAILED'
      WRITE(kuso,*) '  '
      RETURN
      END SUBROUTINE CONFILE
