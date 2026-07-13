SUBROUTINE UPDCAL(KD0,KM0,KY0,KINC,KD1,KM1,KY1,KLMO,KULOUT)
! (C) Copyright 1991- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *UPDCAL*

!     PURPOSE.
!     --------

!     Updates the calendar values.

!**   INTERFACE.
!     ----------

!     CALL UPDCAL(KD0,KM0,KY0,KINC,KD1,KM1,KY1,KLMO,KULOUT)
!          KD0,KM0,KY0 : initial date
!          KINC        : number of days to increment
!          KD1,KM1,KY1 : final date
!          KLMO        : length of the 12 months
!          KULOUT      : output unit (If negative, does not write)

!     METHOD.
!     -------

!     Increases day by day the date, updates if necessary the month and the
!     year.

!     EXTERNALS.
!     ----------

!         NONE

!     AUTHORS.
!     --------

!         M. DEQUE  JAN 91 .
!         Modified 96-07 M. DEQUE Backward increment. of date (Nudging)
!         Modified 99-10 P. VITERBO Control of printing, 2100 compatible
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM) :: KD0,KD1,KINC,KM0,KM1,KULOUT,KY0,KY1

INTEGER(KIND=JPIM) :: KLMO(12)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JD, JM


!*
!     1. LENGTH OF THE MONTHS.
!     ------------------------


DO JM=1,12
  KLMO(JM)=31
  IF(JM == 4.OR.JM == 6.OR.JM == 9.OR.JM == 11)KLMO(JM)=30
  IF(JM == 2)THEN
    IF(MOD(KY0,4) == 0 .AND. MOD(KY0,400) /= 100 &
     &.AND. MOD(KY0,400) /= 200 .AND. MOD(KY0,400) /= 300)THEN
      KLMO(JM)=29
    ELSE
      KLMO(JM)=28
    ENDIF
  ENDIF
ENDDO
KD1=KD0
KM1=KM0
KY1=KY0

!*
!     2. LOOP ON THE DAYS.
!     --------------------


IF(KINC >= 0) THEN
  DO JD=1,KINC
    KD1=KD1+1
    IF(KD1 <= KLMO(KM1)) CYCLE
    KD1=1
    KM1=KM1+1
    IF(KM1 <= 12) CYCLE
    KM1=1
    KY1=KY1+1
    KLMO(2)=28
    IF(MOD(KY1,4) == 0 .AND. MOD(KY1,400) /= 100 &
     &.AND. MOD(KY1,400) /= 200 .AND. MOD(KY1,400) /= 300)KLMO(2)=29
  ENDDO
ELSE
  DO JD=1,-KINC
    KD1=KD1-1
    IF(KD1 > 0) CYCLE
    KM1=1+MOD(KM1+10,12)
    KLMO(2)=28
    IF(MOD(KY1,4) == 0 .AND. MOD(KY1,400) /= 100 &
     &.AND. MOD(KY1,400) /= 200 .AND. MOD(KY1,400) /= 300)KLMO(2)=29
    KD1=KLMO(KM1)
    IF(KM1 == 12)KY1=KY1-1
  ENDDO
ENDIF
IF (KULOUT >=0) WRITE(KULOUT,'("  DATE=",3I5)')KY1,KM1,KD1
RETURN
END SUBROUTINE UPDCAL
