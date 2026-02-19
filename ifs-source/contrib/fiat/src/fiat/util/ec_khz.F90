! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE EC_KHZ(KOREID,KHZ)
USE EC_PARKIND, ONLY : JPIM
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)  :: KOREID
INTEGER(KIND=JPIM),INTENT(OUT) :: KHZ
INTEGER(KIND=JPIM) :: ITHISCORE, ISTAT
INTEGER, external :: EC_COREID ! from ec_env.c
!  /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq
CHARACTER(LEN=80) :: CLSYS
IF (KOREID >= 0) THEN
   ITHISCORE = KOREID
ELSE
   ITHISCORE = EC_COREID()
ENDIF
WRITE(CLSYS,'(A,I0,A)') '/sys/devices/system/cpu/cpu',ITHISCORE,'/cpufreq/scaling_cur_freq'
OPEN(505,FILE=trim(CLSYS),IOSTAT=ISTAT,STATUS='old',ACTION='read')
IF (ISTAT == 0) THEN
   READ(505,*,IOSTAT=ISTAT) KHZ
   CLOSE(505)
ENDIF
IF (ISTAT /= 0) KHZ = 0
END SUBROUTINE EC_KHZ
