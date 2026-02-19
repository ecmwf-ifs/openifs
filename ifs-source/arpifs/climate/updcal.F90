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
SUBROUTINE UPDCAL(KD0,KM0,KY0,KINC,KD1,KM1,KY1,KLMO,KULOUT)

!**** *UPDCAL*

!     PURPOSE.
!     --------

!     Updates the calendar values.

!**   INTERFACE.
!     ----------

!     CALL UPDCAL(...)

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
!      M. DEQUE  JAN 91 .

!     MODIFICATIONS:
!     --------------
! ----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KD0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KM0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KY0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINC 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KD1 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KM1 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KY1 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLMO(12) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM) :: JD, JM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*
!     1. LENGTH OF THE MONTHS.
!     ------------------------

IF (LHOOK) CALL DR_HOOK('UPDCAL',0,ZHOOK_HANDLE)
DO JM=1,12
  KLMO(JM)=31
  IF(JM == 4.OR.JM == 6.OR.JM == 9.OR.JM == 11)KLMO(JM)=30
  IF(JM == 2)THEN
    IF(MOD(KY0,4) == 0 .AND. MOD(KY0,400) /= 100 &
       & .AND. MOD(KY0,400) /= 200 .AND. MOD(KY0,400) /= 300)THEN  
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
     & .AND. MOD(KY1,400) /= 200 .AND. MOD(KY1,400) /= 300)KLMO(2)=29  
  ENDDO
ELSE
  DO JD=1,-KINC
    KD1=KD1-1
    IF(KD1 > 0) CYCLE
    KM1=1+MOD(KM1+10,12)
    KLMO(2)=28
    IF(MOD(KY1,4) == 0 .AND. MOD(KY1,400) /= 100 &
     & .AND. MOD(KY1,400) /= 200 .AND. MOD(KY1,400) /= 300)KLMO(2)=29  
    KD1=KLMO(KM1)
    IF(KM1 == 12)KY1=KY1-1
  ENDDO
ENDIF
IF (KULOUT >=0) WRITE(KULOUT,'("  DATE=",3I5)')KY1,KM1,KD1
IF (LHOOK) CALL DR_HOOK('UPDCAL',1,ZHOOK_HANDLE)
END SUBROUTINE UPDCAL
