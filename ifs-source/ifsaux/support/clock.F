C (C) Copyright 2014- ECMWF.
C
C This software is licensed under the terms of the Apache Licence Version 2.0
C which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
C
C In applying this licence, ECMWF does not waive the privileges and immunities
C granted to it by virtue of its status as an intergovernmental organisation
C nor does it submit to any jurisdiction.

      FUNCTION CLOCK()
!
!     clock() returns the amount of CPU time (in microseconds) used since
!     the first call to clock().  The time reported is the sum of the user
!     and system times of the calling process and its terminated child
!     processes for which it has executed wait() or system() (see wait(2)
!     and system(3S)).  To determine the time in seconds, the value returned
!     by clock() should be divided by the value of the macro CLOCKS_PER_SEC.
!
!     INTRODUCE A PROPER CLCOK ROUTINE HERE
!
!      IMPLICIT LOGICAL (L)
      USE PARKIND1, ONLY : JPIM, JPRB
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: CLOCK
      INTEGER(KIND=JPIM) :: IC, ICR, ICM
      CALL SYSTEM_CLOCK(IC,ICR,ICM)
      IF(ICM.EQ.0) THEN
        PRINT *,' no system clock'
        CLOCK=-99
      ELSE
#ifdef VPP
        CLOCK=IC
#else
        CLOCK=INT(1000.0_JPRB*REAL(IC,KIND=JPRB)/ICR)
#endif
      ENDIF
      RETURN
      ENDFUNCTION CLOCK
