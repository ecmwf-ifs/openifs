! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE SUMPOUT(CDNODE)

!**** *SUMPOUT*   - Additional printings for basic distributed memory environment

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUMPOUT

!        explicit arguments :
!        --------------------
!              CDNODE : output filename prefix. Should depend of nulout to avoid concurrent disk access in write mode

!        implicit arguments :
!        --------------------

!     method.
!     -------
!        see documentation

!     externals.
!     ----------

!     reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     author.
!     -------
!        MPP Group *ECMWF*

!     modifications.
!     --------------
!        original : 95-10-01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib : 26-Feb-2008 Portability fix
!        K. Yessad (Oct 2013): cleanings; calculation of LOUTPUT moved in SUMPINI
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        R. El Khatib : 01-Dec-2014 (re)move call to DATE_AND_TIME
!        R. El Khatib : 18-Mar-2015 EC_LISTING_PATH
!        R. El Khatib : 02-Aug-2018 CDNODE
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : LOUTPUT, MYPROC, MYSETA, MYSETB
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER (LEN = *), INTENT(IN) :: CDNODE

CHARACTER (LEN = 6) :: CLOUT
CHARACTER (LEN = 256) ::  CLENV
LOGICAL :: LL_OPEN
INTEGER(KIND=JPIM) :: IUNIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMPOUT',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF (LOUTPUT) THEN
  CALL GET_ENVIRONMENT_VARIABLE('EC_LISTING_PATH', CLENV)
  WRITE(CLOUT,'(I3.3,A1,I2.2)') MYSETA,'_',MYSETB
  IF (CLENV == ' ') THEN
    OPEN (UNIT=NULOUT, FILE=TRIM(CDNODE)//'.'//CLOUT)  !,ACTION='WRITE',STATUS='REPLACE')
  ELSE
    OPEN (UNIT=NULOUT, FILE=TRIM(CLENV)//'/'//TRIM(CDNODE)//'.'//CLOUT)  !,ACTION='WRITE',STATUS='REPLACE')
  ENDIF
  WRITE(NULOUT,'(/A,I4,A)')&
   & '--- Start of IFS output from processor ',&
   & MYPROC,' ----------------------------------'  
ELSE

! From the Fortran 95 standard, Section 9.3.4:
!
!   If a file is already connected to a unit, execution of an OPEN
!   statement on that file and a different unit is not permitted.
!
! From the Fortran 77 Standard, Section 12.3.2:
!
!   A unit must not be connected to more than one file at the same time,
!   and a file must not be connected to more than one unit at the same time.

  INQUIRE (FILE='/dev/null', OPENED=LL_OPEN, NUMBER=IUNIT)
  IF (LL_OPEN) CLOSE(IUNIT)
  OPEN(UNIT=NULOUT, FILE='/dev/null')
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMPOUT',1,ZHOOK_HANDLE)
END SUBROUTINE SUMPOUT
