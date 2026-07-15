! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WPOSNAM(KULNAM, CDNAML, LDEOF)

!**** *WPOSNAM* - position namelist file for reading

!     Purpose.
!     --------
!     The sole purpose of this is to fix a bug of the SGI f90
!     compiler, which does not allow to trap an EOF condition
!     while reading a namelist.
!     So we do it here with an ordinary read.

!**   Interface.
!     ----------
!        *CALL* *WPOSNAM*(..)

!        Explicit arguments :     KULNAM - file unit number (input)
!        --------------------     CDNAML - namelist name    (input)
!                                 LPEOF  - end of file switch (output)


!        Implicit arguments :     None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Bjorn Hansen November 1998  *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-06-22 from Mats Hamrud as posnam
!     --------------------------------------------------------------

      USE PARKIND_WAVE,  ONLY : JWIM
      USE YOWABORT, ONLY : WAM_ABORT
      USE EC_LUN   , ONLY : NULERR

!     --------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL, INTENT(OUT)      :: LDEOF
      INTEGER(KIND=JWIM), INTENT(IN)       :: KULNAM
      CHARACTER(LEN=*), INTENT(IN) :: CDNAML

      INTEGER(KIND=JWIM) :: ILEN, IND1

! === END OF INTERFACE BLOCK ===

      CHARACTER(LEN=120) :: CLINE
      CHARACTER(LEN=1) :: CLTEST

!*       1.    POSITION FILE
!              -------------
      LDEOF = .FALSE.

      ILEN=LEN(CDNAML)
 102  CONTINUE
      CLINE=' '
      READ(KULNAM,'(A)', END=9000, ERR=9001) CLINE
      IND1=INDEX(CLINE,'&'//CDNAML)
      IF(IND1.EQ.0) GO TO 102
      CLTEST=CLINE(IND1+ILEN+1:IND1+ILEN+1)
      IF((LGE(CLTEST,'0').AND.LLE(CLTEST,'9')).OR.                      &
     &   (LGE(CLTEST,'A').AND.LLE(CLTEST,'Z'))) GO TO 102
      BACKSPACE(KULNAM)


!     ------------------------------------------------------------------

      RETURN

!     ------------------------------------------------------------------

 9000 CONTINUE

      LDEOF = .TRUE.
      RETURN

 9001 CONTINUE

      WRITE(NULERR,*) '*********************'
      WRITE(NULERR,*) 'READ ERROR IN WPOSNAM ' 
      WRITE(NULERR,*) '*********************'
      CALL WAM_ABORT("READ ERROR IN WPOSNAM",__FILENAME__,__LINE__)

      RETURN

      END SUBROUTINE WPOSNAM
