! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GSTATS_PSUT

! MEASURE PARALLELL START UP TIME

USE EC_PARKIND  ,ONLY : JPRD, JPIM

USE YOMGSTATS
USE MPL_MODULE


IMPLICIT NONE

INTEGER(KIND=JPIM) :: ITAG,ILSEND,ILBUF,JROC,IRECV,ISEND
REAL(KIND=JPRD) :: ZTBUF(2),ZCLOCK,ZCLOCKB

#include "user_clock.intfb.h"

IF(NPROC_STATS > 1) THEN
  CALL USER_CLOCK(PELAPSED_TIME=ZCLOCKB)
  CALL MPL_BARRIER
  CALL USER_CLOCK(PELAPSED_TIME=ZCLOCK)
  ITAG = JPTAGSTAT
  IF (MYPROC_STATS == 1 ) THEN
    ALLOCATE(TIME_START(NPROC_STATS))
    TIME_START(1) = ZCLOCKB - TIMELCALL(0)
    ILBUF = 2
  ENDIF
  DO JROC=2,NPROC_STATS
   IF (MYPROC_STATS .eq. jroc ) THEN
    ZTBUF(1) = ZCLOCKB
    ZTBUF(2) = ZCLOCK
    ILSEND = 2
    ISEND = 1
    CALL MPL_SEND(ZTBUF(1:ILSEND),KDEST=NPRCIDS_STATS(ISEND), &
     & KTAG=ITAG,CDSTRING='SUSTATS:')
   ELSEIF (MYPROC_STATS == 1 ) THEN
      IRECV = JROC
      CALL MPL_RECV(ZTBUF(1:ILBUF),KSOURCE=NPRCIDS_STATS(IRECV), &
       & KTAG=ITAG,CDSTRING='SUSTATS:')
      TIME_START(JROC) = ZTBUF(1) - TIMELCALL(0) -(ZTBUF(2)-ZCLOCK)
   ENDIF
   CALL MPL_BARRIER
  ENDDO
ENDIF

END SUBROUTINE GSTATS_PSUT
