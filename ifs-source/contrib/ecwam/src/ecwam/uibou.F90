! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE UIBOU (IU05, IU06, CDATEA, CDATEE, IDELPRF,            &
     &                  CDATES, IDELFI, USERID,                         &
     &                  RUNDI, FILEDI, PATHI, RUNDO, FILEDO, PATHO,     &
     &                  LREAL, LSWAN)

! -----------------------------------------------------------------
 
!**** *UIBOU* - ROUTINE TO READ USER INPUT FOR PROG BOUINT.
 
!     H. GUNTHER     GKSS/ECMWF     JANUARY 1991

!!!! it should be replaced by a proper namelist !!!!
 
!*    PURPOSE.
!     --------
 
!       READ USER INPUT OF PROGRAM BOUINT
 
!**   INTERFACE.
!     ----------
 
!       *CALL* *UIBOU (IU05, IU06, IDATEA, CDATEE, IDELPRF,
!                      CDATES, IDELFI, USERID,
!                      RUNDI, FILEDI, PATHI, RUNDO, FIEDO, PATHO,
!                      LREAL, LSWAN)
!         *IU05*    INTEGER    INPUT UNIT.
!         *IU06*    INTEGER    PRINTER OUTPUT UNIT.
!         *CDATEA*  CHARACTER  START DATE OF OUTPUT (YYYYMMDDHHMMSS).
!         *CDATEE*  INTEGER    END DATE OF OUTPUT (YYYYMMDDHHMMSS).
!         *IDELPRF* INTEGER    OUTPUT TIME INCREMENT (SECONDS).
!         *USERID*  CHARACTER  USERID OF INPUT FILES.
!         *IDATES*  INTEGER    LAST FIELD DATE IN FIRST INPUT FILE
!                              (YYMMDDHHMM).
!         *IDELFI*  INTEGER    INPUT FILE DATE INCREMENT (SECONDS).
!         *USERID*  CHARACTER  USERID OF INPUT FILES.
!         *RUNDI*   CHARACTER  RUN  ID OF INPUT FILES.
!         *FILEDI*  CHARACTER  FILE ID OF INPUT FILES.
!         *PATHI *  CHARACTER  DIRECTORY OF INPUT FILES.
!         *RUNDO*   CHARACTER  RUN  ID OF OUTPUT FILES.
!         *FILEDO*  CHARACTER  FILE ID OF OUTPUT FILES.
!         *PATHO *  CHARACTER  DIRECTORY OF OUTPUT FILES.
!         *LREAL*   LOGICAL    IF TRUE REAL*4 OUTPUT.
!         *LSWAN*   LOGICAL    IF TRUE SWAn LIKE OUTPUT.
 
!     METHOD.
!     -------
 
!        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
!         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
!         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
!            (RIGHT-JUSTIFIED, BUT CHARACTER LEFT-JUSTIFIED)
 
!     EXTERNALS.
!     ----------
 
!       *ABORT*     - TERMINATES PROCESSING.
 
!     REFERENCE.
!     ----------
 
!       NONE.
 
! ------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ------------------------------------------------------------------

      IMPLICIT NONE 
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: IU05, IU06, IDELPRF, IDELFI  
      INTEGER(KIND=JWIM) :: ICOUNT, IOS
 
      CHARACTER(LEN=3) :: USERID, RUNDI, FILEDI, RUNDO, FILEDO 
      CHARACTER(LEN=14) :: CDATEA, CDATEE, CDATES
      CHARACTER(LEN=60) :: PATHI, PATHO  
      CHARACTER(LEN=120) :: LINE

      LOGICAL :: LREAL, LSWAN
 
! --------------------------------------------------------------------
 
!*    1. READ USER INPUT
!        ---------------
 
      LREAL=.FALSE.
      LSWAN=.FALSE.

      ICOUNT = 0
 1000 CONTINUE
      READ(IU05, '(A120)',ERR=4000,END=4000,IOSTAT=IOS) LINE
      IF (LINE(1:1).EQ.'C') GOTO 1000
      ICOUNT = ICOUNT + 1
      IF (ICOUNT.EQ.1) THEN
         CDATEA=LINE( 2:15)
         CDATEE=LINE(18:31)
         READ(LINE(34:40),'(I7)',ERR=4100,IOSTAT=IOS) IDELPRF
         IF (LINE(42:42).EQ.'H') IDELPRF = IDELPRF*3600
         GOTO 1000
      ELSE IF (ICOUNT.EQ.2) THEN
!         READ(LINE( 2:11),'(I10)',ERR=4100,IOSTAT=IOS) IDATES
         CDATES=LINE( 2:15)
         READ(LINE(18:24),'(I7)',ERR=4100,IOSTAT=IOS) IDELFI
         IF (LINE(26:26).EQ.'H') IDELFI = IDELFI*3600
         USERID = LINE(29:31)
         GOTO 1000
      ELSE IF (ICOUNT.EQ.3) THEN
         RUNDI  = LINE(2:4)
         FILEDI = LINE(7:9)
         PATHI  = LINE(12:71)
         GOTO 1000
      ELSE IF (ICOUNT.EQ.4) THEN
         RUNDO  = LINE(2:4)
         FILEDO = LINE(7:9)
         PATHO  = LINE(12:71)
         GOTO 1000
      ELSE IF (ICOUNT.EQ.5) THEN
         READ(LINE(2:2),'(L1)') LREAL
         READ(LINE(4:4),'(L1)') LSWAN
      ENDIF
 
!     2. PRINT USER INPUT
!        ----------------
 2000 CONTINUE
      WRITE(IU06,*) ' USER INPUT:'
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' START  DATE (FORMAT:YYMMDDHHMM) : ',CDATEA,       &
     &           ' END DATE :',CDATEE
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' OUTPUT EVERY ',IDELPRF ,' SECONDS'
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' INPUT DATA WILL BE FETCHED WITH:'
      WRITE(IU06,*) ' USER ID IS ..................... ', USERID
      WRITE(IU06,*) ' RUN  ID IS ..................... ', RUNDI
      WRITE(IU06,*) ' FILE ID IS ..................... ', FILEDI
      WRITE(IU06,*) ' DIRECTORY NAME IS .............. ', PATHI
      WRITE(IU06,*) ' DATE OF LAST FIELD IN FIRST FILE ', CDATES
      WRITE(IU06,*) ' A NEW FILE WILL BE FETCHED EVERY ', IDELFI,       &
     &              ' SECONDS'
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' OUTPUT DATA WILL BE STORED WITH:'
      WRITE(IU06,*) ' USER ID IS ..................... ', USERID
      WRITE(IU06,*) ' RUN  ID IS ..................... ', RUNDO
      WRITE(IU06,*) ' FILE ID IS ..................... ', FILEDO
      WRITE(IU06,*) ' DIRECTORY NAME IS .............. ', PATHO
      WRITE(IU06,*) ' DATE OF LAST FIELD IN FIRST FILE ', CDATES
      WRITE(IU06,*) ' A NEW FILE WILL BE FETCHED EVERY ', IDELFI,       &
     &              ' SECONDS'
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' LREAL = ', LREAL
      WRITE(IU06,*) ' LSWAN = ', LSWAN
 
!*    3. CHECK CONSISTENCY OF INPUT DATA
!        -------------------------------
 
      IF (CDATEE.LT.CDATEA) THEN
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR IN SUB. UIBOU            *'
         WRITE(IU06,*) '*    =========================            *'
         WRITE(IU06,*) '* END DATE IS BEFORE START DATE           *'
         WRITE(IU06,*) '* START DATE = ', CDATEA
         WRITE(IU06,*) '* END  DATE  = ', CDATEE
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* CORRECT USER INPUT                      *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '* ---------------   --------------        *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
      ENDIF

      RETURN
 
!*    4. READ ERROR MESSAGES
!        -------------------
 
 4000 CONTINUE
         WRITE(IU06,*) ' ********************************************'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *     FATAL ERROR IN SUB. UIBOU            *'
         WRITE(IU06,*) ' *     =========================            *'
         WRITE(IU06,*) ' * READ ERROR ON INPUT FILE:                *'
         WRITE(IU06,*) ' * ERROR IS LATER THAN ICOUNT = ', ICOUNT
         WRITE(IU06,*) ' * LAST LINE READ IS     LINE = ', LINE
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS         *'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' ********************************************'
         CALL ABORT1
 4100 CONTINUE
         WRITE(IU06,*) ' ********************************************'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *     FATAL ERROR IN SUB. UIBOU            *'
         WRITE(IU06,*) ' *     =========================            *'
         WRITE(IU06,*) ' * READ ERROR ON CHARACTER STRING           *'
         WRITE(IU06,*) ' * ERROR IS IN DATA LINE ICOUNT = ', ICOUNT
         WRITE(IU06,*) ' * CHARACTER STRING IS   LINE = ', LINE
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS         *'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' ********************************************'
         CALL ABORT1

      END SUBROUTINE UIBOU
