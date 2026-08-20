! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WGRIBOUT (IU06, ITEST, LFDB, IU, IGRIB_HANDLE, ISIZE, KGRIB_BUFR)
! ----------------------------------------------------------------------

!****  *WGRIBOUT* OUTPUT GRIB MESSAGE. 

!       PURPOSE.
!       --------

!      OUTPUT GRIB MESSAGE TO FDB (ECMWF) OR TO FILE.

!**    INTERFACE.
!      ----------
!        *CALL* *WGRIBOUT (IU06, ITEST, LFDB, IU, IGRIB_HANDLE, ISIZE, KGRIB_BUFR)
!          *IU06*    LOGFILE OUTPUT UNIT.
!          *ITEST*   TEST OUTPUT GIVEN IF ITEST GT 2.
!          *LFDB*    LOGICAL SWITCH TO ACTIVATE FDB PROCESSING
!                    .T. TO ACTIVATE, .F. TO IGNORE
!          *IU*      UNIT TO WRITE O FILE
!                    !!! IT WILL NEED TO BE OPEN AN
!                    !!! WITH GRIB_API SOFTWARE !!!
!          *IGRIB_HANDLE* GRIB MESSAGE HANDLE.
!          *ISIZE*   GRIB MESSAGE SIZE.
!          *KGRIB_BUFR* GRIB MESSAGE. 

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIB  , ONLY : JPKSIZE_T, &
                          & IGRIB_GET_VALUE, &
                          & IGRIB_GET_MESSAGE_SIZE, &
                          & IGRIB_WRITE_BYTES
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wgrib2fdb.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IGRIB_HANDLE, ISIZE
      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ITEST, IU
      INTEGER(KIND=JWIM), DIMENSION(ISIZE), INTENT(INOUT) :: KGRIB_BUFR
      LOGICAL, INTENT(IN) :: LFDB

      INTEGER :: IERR, ITABPAR, ICLASS
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=4) :: CEXPVER
      CHARACTER(LEN=12) :: C12


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WGRIBOUT',0,ZHOOK_HANDLE)

      CALL GSTATS(1709,0)

!*    1. SAVE ENCODED RESULT
!        -------------------

      IF(ITEST.GT.0) THEN
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'paramId',ITABPAR)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'expver',C12)
        CEXPVER=C12(1:4)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',C12)
        CSTREAM=C12(1:4)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'type',ICLASS)
        WRITE(IU06,*)'  '
        WRITE(IU06,*)'   SUB. WGRIBOUT : PARAMETER ',ITABPAR,           &
     &                                 ' EXPVER=',CEXPVER,              &
     &                                 ' STREAM= ', CSTREAM,            &
     &                                 ' CLASS = ', ICLASS
        CALL FLUSH(IU06)
      ENDIF


      IF(LFDB) THEN
!*     2.1 WRITE DATA TO FDB
!          -----------------
        CALL WGRIB2FDB (IU06, ITEST,                                    &
     &                  IGRIB_HANDLE, ISIZE, KGRIB_BUFR,                &
     &                  IERR)
        IF(ITEST.GT.0) THEN
          WRITE(IU06,*) '   SUB. WGRIBOUT: GRIB DATA WRITTEN TO FDB'
          CALL FLUSH(IU06)
        ENDIF

        IF(IERR.NE.0)THEN
          WRITE(IU06,*) ' ------------------------'
          WRITE(IU06,*) ' ERROR ACCESSING FDB '
          WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
          WRITE(IU06,*) ' ------------------------'
          WRITE(NULERR,*) ' ------------------------'
          WRITE(NULERR,*) ' ERROR ACCESSING FDB '
          WRITE(NULERR,*) ' FDB ERROR CODE IS ',IERR
          WRITE(NULERR,*) ' ------------------------'
          CALL ABORT1
        ENDIF

      ELSE

!*      2.2 WRITE PACKED DATA TO FILE CONNECTED TO UNIT IU.
!           ----------------------------------------------

        CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
        CALL IGRIB_WRITE_BYTES(IU,KGRIB_BUFR,KBYTES)

        IF (ITEST.GT.2) THEN
            WRITE(IU06,*) '   SUB. WGRIBOUT : DATA WRITTEN TO FILE'
        ENDIF
      ENDIF

      CALL GSTATS(1709,1)

      IF (LHOOK) CALL DR_HOOK('WGRIBOUT',1,ZHOOK_HANDLE)

      END SUBROUTINE WGRIBOUT
