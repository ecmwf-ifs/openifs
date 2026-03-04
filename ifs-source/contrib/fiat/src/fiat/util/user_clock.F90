! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE USER_CLOCK(PELAPSED_TIME,PELAPSED_TIME_SINCE,PVECTOR_CP,PTOTAL_CP)

!**** *USER_CLOCK* - interface to system dependent timer routines

!     Purpose.
!     --------
!        Returns elapsed and CP from the start of execution.
!        Elapsed time is made relative to the first call to USER_CLOCK.

!**   Interface.
!     ----------
!        ZTIME=USER_CLOCK(PELAPSED_TIME,PELAPSED_TIME_SINCE,
!                         PVECTOR_CP,PTOTAL_CP)

!        Explicit arguments: (All are optional arguments)
!                           PELAPSED_TIME=wall clock time (seconds)
!                           PELAPSED_TIME_SINCE=wall clock time (seconds)
!                             change from input value of this parameter
!                           PVECTOR_CP=CP vector time  (seconds)
!                           PTOTAL_CP=total CP time   (seconds)

!     Author.
!     -------
!        D.Dent      *ECMWF*

!     External References:
!     -------------------

!        TIMEF,CPTIME

!     Modifications.
!     --------------
!        Original  : 97-09-25
!      F. Vana  05-Mar-2015  Support for single precision
!     ----------------------------------------------------------


USE EC_PARKIND  ,ONLY : JPRD, JPIM

IMPLICIT NONE

REAL(KIND=JPRD),INTENT(OUT) :: PELAPSED_TIME,PVECTOR_CP,PTOTAL_CP
REAL(KIND=JPRD),INTENT(INOUT) :: PELAPSED_TIME_SINCE
OPTIONAL            PELAPSED_TIME,PELAPSED_TIME_SINCE
OPTIONAL            PVECTOR_CP,PTOTAL_CP
REAL(KIND=JPRD)      :: ZVECTOR_CP,ZTOTAL_CP,ZWALL
REAL(KIND=JPRD),EXTERNAL :: TIMEF


! === END OF INTERFACE BLOCK ===
IF(PRESENT(PELAPSED_TIME).OR. PRESENT(PELAPSED_TIME_SINCE)) THEN

  ZWALL=TIMEF()
!             TIMEF returns milliseconds since first call to TIMEF
  IF(PRESENT(PELAPSED_TIME)) THEN
    PELAPSED_TIME=ZWALL*1.0E-3_JPRD
  ENDIF
  IF(PRESENT(PELAPSED_TIME_SINCE)) THEN
    PELAPSED_TIME_SINCE=ZWALL*1.0E-3_JPRD - PELAPSED_TIME_SINCE
  ENDIF
ENDIF

IF( PRESENT(PVECTOR_CP) .OR. PRESENT(PTOTAL_CP) ) THEN
  CALL CPTIME(ZVECTOR_CP,ZTOTAL_CP)
ENDIF
IF( PRESENT(PVECTOR_CP) ) THEN
  PVECTOR_CP=ZVECTOR_CP
ENDIF
IF( PRESENT(PTOTAL_CP) ) THEN
  PTOTAL_CP=ZTOTAL_CP
ENDIF

RETURN
END SUBROUTINE USER_CLOCK



