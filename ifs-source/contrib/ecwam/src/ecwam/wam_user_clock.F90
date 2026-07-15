! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      FUNCTION WAM_USER_CLOCK()

!   Returns system clock converted to micro-seconds

      USE PARKIND_WAVE, ONLY : JWIM, JWRU

      IMPLICIT NONE
      REAL(JWRU) :: WAM_USER_CLOCK
      INTEGER(KIND=JWIM) :: ICOUNT,ICOUNT_RATE,ICOUNT_MAX

      CALL SYSTEM_CLOCK(ICOUNT,ICOUNT_RATE,ICOUNT_MAX)
      IF(ICOUNT.NE.HUGE(0)) THEN
        WAM_USER_CLOCK = REAL(ICOUNT,JWRU) * 1.E6_JWRU /            &
     &                   REAL(ICOUNT_RATE,JWRU)
      ELSE
        WAM_USER_CLOCK = 0.0_JWRU
      ENDIF    

      RETURN
      END FUNCTION WAM_USER_CLOCK




