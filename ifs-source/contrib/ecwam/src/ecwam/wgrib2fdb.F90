! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!=======================================================================
      SUBROUTINE WGRIB2FDB (KUSO, KTEST,                                &
     &                      IGRIB_HANDLE, KLEN, KGRIB,                  &
     &                      KERR)
!=======================================================================
!**
!**** NAME    *WGRIB2FDB*
!**** ----
!**
!**   PURPOSE
!**   -------
!**      WRITE FIELD TO THE FIELD DATA BASE.
!**
!**   INTERFACE
!**   ---------
!**      INPUT:*KUSO*         LOGICAL UNIT FOR STANDARD OUTPUT.
!**            *KTEST*        SWITCH DIAGNOSTICS OUTPUT ON IF KTEST GT 1.
!**            *IGRIB_HANDLE* GRIB HANDLE CONTAINING THE ENCODED DATA.
!**            *KLEN*         GRIB MESSAGE LENGTH.
!**            *KGRIB*        GRIB MESSAGE.
!**            *KERR*         0 IF NO ERRORS ENCOUNTERED.
!**
!**   METHOD
!**   ------
!**      FIELD DATA BASE ROUTINES ARE USED TO OPEN THE FIELD
!**      DATA BASE AND FILE, AND TO WRITE/READ THE FIELD INTO/
!**      FROM THE FIELD DATA BASE.
!**
!**
!**   AUTHOR
!**   ------
!**     JEAN BIDLOT   ECMWF   AUGUST 2009. 

!     ------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPRECI

      USE WAM_MULTIO_MOD, ONLY : WAM_MULTIO_WRITE, WAM_MULTIO_FLUSH
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------
      IMPLICIT NONE
#include "wam_u2l1cr.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KUSO, KTEST, IGRIB_HANDLE
      INTEGER(KIND=JWIM), INTENT(IN) :: KLEN 
      INTEGER(KIND=JWIM), INTENT(OUT) :: KERR
      INTEGER(KIND=JWIM), DIMENSION(KLEN), INTENT(INOUT) :: KGRIB

      INTEGER(KIND=JWIM) :: ISTAT

!!      INTEGER(KIND=JPKSIZE_T) :: KBYTES
!!      INTEGER(KIND=JWIM), SAVE :: IFILE_HANDLE = -999
!!      CHARACTER(LEN=64) :: CLFILE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WGRIB2FDB',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!for debugging
!!      CLFILE = 'allfields_%p.grib'
!!      CALL EXPAND_STRING(IRANK,NPROC,0,0,CLFILE,1)
!!      IF(IFILE_HANDLE == -999) THEN
!!        CALL IGRIB_OPEN_FILE(IFILE_HANDLE,TRIM(CLFILE),'w')
!!      ENDIF


!     THIS VALUE WILL BE RETURNED IF EVERYTHING GOES OK.
      KERR = 0

!     ------------------------------------------------------------------

!*    3.  GET AND WRITE FIELD TO FDB

!     ------------------------------------------------------------------

!for debugging
!!      write(*,*) 'also writting to file the data destined to fdb !!!'
!!      CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
!!      CALL IGRIB_WRITE_BYTES(IFILE_HANDLE,KGRIB,KBYTES)

      CALL WAM_MULTIO_WRITE( KGRIB, KLEN, ISTAT )
      IF ( ISTAT .NE. 0 ) THEN
        WRITE(KUSO, '( "Error\ /WGRIB2FDB/ Failed to write to fdb")' )
        CALL WAM_MULTIO_FLUSH()
        KERR = 1
        RETURN
      ENDIF
      IF (KTEST.GT.1) WRITE(KUSO,'("\ /WGRIB2FDB/ iwritefdb klen:", i10, &
     &                      " status:", i4)') KLEN, ISTAT

      IF (LHOOK) CALL DR_HOOK('WGRIB2FDB',1,ZHOOK_HANDLE)

      END SUBROUTINE WGRIB2FDB
