! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      REAL FUNCTION WAM_USER_CLOCK()

!   Returns system clock converted to micro-seconds

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: ICOUNT,ICOUNT_RATE,ICOUNT_MAX

      CALL SYSTEM_CLOCK(ICOUNT,ICOUNT_RATE,ICOUNT_MAX)
      IF(ICOUNT.NE.HUGE(0)) THEN
        WAM_USER_CLOCK = REAL(ICOUNT+0.0_JWRB) * 1.E6_JWRB /            &
     &                   REAL(ICOUNT_RATE+0.0_JWRB)
      ELSE
        WAM_USER_CLOCK = 0.0_JWRB
      ENDIF    

      RETURN
      END FUNCTION WAM_USER_CLOCK




