! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ----------------------------------------------------------------------

      SUBROUTINE FILE_TRANSFER (IU06, IUNIT, FILENA, CDOPT)

! ----------------------------------------------------------------------

!**** *FILE_TRANSFER* - MOVE FORTRAN UNIT FILE TO A FILE 

!     THIS IS QUITE OLD FASHIONED. ALL WRITES TO FILE SHOULD BE
!     REWRITTEN !!! AT ECMWF, WE NORMALLY DO NOT USE THIS OPTION
!     AS WE WRITE IN GRIB FORMAT DIRECTLY INTO THE FIELD DATA
!     BASE !!! (SORRY).

!     J. BIDLOT     ECMWF JUNE  1996  MESSAGE PASSING
!     B. HANSEN     ECMWF APRIL 1996  MOVE CYMPPATH TO COMTEXT AND ADD 
!                                     CDOPT.

!*    PURPOSE.
!     --------

!     MOVE FORTRAN UNIT FILE TO A FILE. 

!**   INTERFACE.
!     ----------

!       *CALL* *FILE_TRANSFER*(IU06, IUNIT, FILENA) 

!         *IU06*   INTEGER    UNIT FOR PRINTER MESSAGES.
!         *IUNIT*  INTEGER    FORTRAN UNIT.
!         *FILENA* CHARACTER  DESTINATION FILENAME
!         *CDOPT*  CHARACTER  OPTION
!                             = 'S' SAVE FILE CONNECTED TO IUNIT
!                             = 'G' GET FILE AND CONNECT FILE TO IUNIT

!     EXTERNALS.
!     ----------

!     METHOD.
!     -------

!     USE A SYSTEM CALL TO MOVE THE FILE. 

!     REFERENCES.
!      -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH


      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: IU06, IUNIT
      INTEGER(KIND=JWIM) :: LNAME, LNF, LNT

      CHARACTER(LEN=1) :: MODE, CDOPT 
      CHARACTER(LEN=*) :: FILENA
      CHARACTER(LEN=160) :: FILENAM, FILENAMEF,  FILENAMET

      LOGICAL :: LLEXIST

! ----------------------------------------------------------------------
      LLEXIST=.FALSE.

      IF (IUNIT.LT.10) THEN
        WRITE(FILENAMEF,'(A,I1)') 'fort.',IUNIT
      ELSEIF (IUNIT.LT.100) THEN
        WRITE(FILENAMEF,'(A,I2)') 'fort.',IUNIT
      ELSE
        WRITE(FILENAMEF,'(A,I3)') 'fort.',IUNIT
      ENDIF
      FILENAM=FILENAMEF
      LNAME = LEN_TRIM(FILENAM)
      INQUIRE(FILE=FILENAM(1:LNAME),EXIST=LLEXIST)
      IF(.NOT. LLEXIST) THEN
        IF(ICPLEN.GT.0 ) THEN
          FILENAMEF=CPATH(1:ICPLEN)//'/'//FILENAM(1:LNAME)
        ENDIF
      ENDIF

      LNF = LEN_TRIM(FILENAMEF)
      LLEXIST=.FALSE.
      INQUIRE(FILE=FILENAMEF(1:LNF),EXIST=LLEXIST)
      IF(.NOT.LLEXIST .AND. CDOPT.NE.'G') THEN
        WRITE (IU06,*) '******************************************'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO INQUIRE      *'
        WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
        WRITE (IU06,*) '*   COULD NOT FIND FILE ',FILENAMEF
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '******************************************'
        CALL ABORT1
      ENDIF

      FILENAM=''
      LNAME = LEN_TRIM(FILENA)
      FILENAM(1:LNAME)=FILENA(1:LNAME)

      LNAME = LEN_TRIM(FILENAM)
      IF(ICPLEN.EQ.0) THEN
        FILENAMET=FILENAM(1:LNAME)
      ELSE
        FILENAMET=CPATH(1:ICPLEN)//'/'//FILENAM(1:LNAME)
      ENDIF

      IF (CDOPT.EQ.'G') THEN
        FILENAM = FILENAMET
        FILENAMET = FILENAMEF
        FILENAMEF = FILENAM
      ENDIF

      WRITE(IU06,*) 'FILE_TRANSFER: MOVE FILE ',                        &
     &              FILENAMEF(1:LEN_TRIM(FILENAMEF)), ' TO ',           &
     &              FILENAMET(1:LEN_TRIM(FILENAMET))

      LNF = LEN_TRIM(FILENAMEF)
      LNT = LEN_TRIM(FILENAMET)
      CALL SYSTEM('mv '//FILENAMEF(1:LNF)//' '//FILENAMET(1:LNT))

      END SUBROUTINE FILE_TRANSFER
