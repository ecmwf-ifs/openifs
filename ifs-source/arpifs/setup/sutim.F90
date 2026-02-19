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

SUBROUTINE SUTIM(KULOUT)

!**** *SUTIM * - Initialize the common containing real time

!     Purpose.
!     --------
!        Initialisation of the common YOMTIM

!**   Interface.
!     ----------
!        *CALL* *SUTIM(KULOUT)*

!        Explicit arguments :
!        --------------------

!           KULOUT : logical unit for the output

!        Implicit arguments :
!        --------------------
!           common YOMTIM

!     Method.
!     -------
!        Use of the facilities of the CRAY timing functions

!     Externals.
!     ----------
!        USER_CLOCK,DATE_AND_TIME

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!   R. El Khatib : 01-Dec-2014 EC_DATE_AND_TIME
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMTIM   , ONLY : RSTART   ,RVSTART  ,RTIMEF

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

INTEGER(KIND=JPIM) :: J

#include "user_clock.intfb.h"

INTEGER(KIND=JPIM) :: IVALUES(8)
CHARACTER (LEN = 10) :: CLDATEOD,CLTIMEOD,CLZONEOD
CHARACTER (LEN = 1) ::  CLENV
LOGICAL :: LLTIMER, LLDATER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1. INITIALIZE CLOCKS.
!           ------------------

IF (LHOOK) CALL DR_HOOK('SUTIM',0,ZHOOK_HANDLE)

CALL GET_ENVIRONMENT_VARIABLE('EC_DATE_AND_TIME', CLENV)
! EC_DATE_AND_TIME =  0 : no date nor time
! EC_DATE_AND_TIME =  1 : time only
! EC_DATE_AND_TIME = -1 : date only
! EC_DATE_AND_TIME =  any other value : date and time
LLTIMER = (CLENV /= '0').AND.(CLENV /= '-1')
LLDATER = (CLENV /= '0').AND.(CLENV /= '1')

IF (LLTIMER) THEN
  CALL DATE_AND_TIME(CLDATEOD,CLTIMEOD,CLZONEOD,IVALUES)
  WRITE(KULOUT,'(A,A,'':'',A,'':'',A/)') ' TIME OF START= ',&
   & CLTIMEOD(1:2),CLTIMEOD(3:4),CLTIMEOD(5:10)
ENDIF

CALL USER_CLOCK(PELAPSED_TIME=RTIMEF,PVECTOR_CP=RVSTART,PTOTAL_CP=RSTART)

IF (LLDATER) THEN
  WRITE(UNIT=KULOUT,FMT='('' ***   Real world time    ***'')')
  WRITE(UNIT=KULOUT,&
   & FMT='('' Date : '',I4.4,''-'',I2.2,''-'',I2.2,'' Time : ''&
   & ,I2.2,'':'',I2.2,'':'',I2.2)')&
   & (IVALUES(J),J=1,3),(IVALUES(J),J=5,7)  
ENDIF
IF (LLTIMER) THEN
  WRITE(KULOUT,'(1X,F17.10,'' s since start of run'')')RSTART
  WRITE(KULOUT,'(1X,F17.10,'' s from last call  '')')RTIMEF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUTIM',1,ZHOOK_HANDLE)
END SUBROUTINE SUTIM
