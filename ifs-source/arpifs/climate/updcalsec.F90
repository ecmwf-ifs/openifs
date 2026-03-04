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

SUBROUTINE UPDCALSEC(KH0,KD0,KM0,KY0,KINC,KINCS,KH1,KMN1,KS1,KD1,KM1,KY1,KLMO,KULOUT)

!**** *UPDCALSEC*

!     PURPOSE.
!     --------

!     Updates the calendar values. In this version KINC and KINCS are positive
!     they describe the lapse time since the initial date/time of the model run 

!**   INTERFACE.
!     ----------

!     CALL UPDCALSEC(...)

!          KD0,KM0,KY0 : initial date
!          KH0         : initial time in hours
!          KINC        : number of days to increment
!          KINCS       : number of seconds to increment (modulo number of days)
!          KD1,KM1,KY1 : final date
!          KH1,KMN1,KS1: final time (hours,minutes,secondes)
!          KLMO        : length of the 12 months
!          KULOUT      : output unit (If negative, does not write)

!     METHOD.
!     -------

!     calculate the new date and the new time using the given increment, updates if necessary the month and the
!     year.

!     EXTERNALS.
!     ----------

!         NONE

!     AUTHORS.
!     --------
!      S. Serrar
! ----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KD0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KM0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KY0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KH0
INTEGER(KIND=JPIM),INTENT(IN)    :: KINC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINCS
INTEGER(KIND=JPIM),INTENT(OUT)   :: KD1 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KM1 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KY1
INTEGER(KIND=JPIM),INTENT(OUT)   :: KH1
INTEGER(KIND=JPIM),INTENT(OUT)   :: KMN1
INTEGER(KIND=JPIM),INTENT(OUT)   :: KS1
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLMO(12) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM) :: JD, JM , INCHOUR, INCMIN 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*
!     1. LENGTH OF THE MONTHS.
!     ------------------------

IF (LHOOK) CALL DR_HOOK('UPDCALSEC',0,ZHOOK_HANDLE)
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
KH1=KH0
KMN1=0_JPIM
KS1=0_JPIM
KD1=KD0
KM1=KM0
KY1=KY0

!*
!     2. LOOP ON THE DAYS.
!     --------------------
! Here KINC is positive !!!!

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

!*
!     3. CALCULATE FINAL TIME
!     -----------------------

! Here KINCS is positive

INCHOUR=INT(KINCS/3600._JPRB)
INCMIN=(KINCS/3600._JPRB-INCHOUR)*60_JPIM

KH1=KH1+INCHOUR
KMN1=KMN1+INCMIN
KS1=KS1+KINCS-INCHOUR*3600_JPIM-INCMIN*60_JPIM 

IF(KH1 >= 24) THEN
  KH1=KH1-24
  KD1=KD1+1
  IF(KD1 > KLMO(KM1)) THEN
     KD1=1
     KM1=KM1+1
  ENDIF
  IF(KM1 > 12) THEN
    KM1=1
    KY1=KY1+1
  ENDIF
ENDIF

IF (KULOUT >=0) WRITE(KULOUT,'("  DATE=",3I5)')KY1,KM1,KD1
IF (KULOUT >=0) WRITE(KULOUT,'("  TIME=",3I5)')KH1,KMN1,KS1

IF (LHOOK) CALL DR_HOOK('UPDCALSEC',1,ZHOOK_HANDLE)
END SUBROUTINE UPDCALSEC
