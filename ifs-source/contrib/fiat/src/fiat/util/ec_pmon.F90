! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE EC_PMON(ENERGY,POWER)
USE EC_PARKIND, ONLY : JPIM, JPIB
IMPLICIT NONE
INTEGER(KIND=JPIB),INTENT(OUT) :: ENERGY,POWER
INTEGER(KIND=JPIB),SAVE :: ENERGY_START = 0
INTEGER(KIND=JPIM),SAVE :: MONINIT = 0
INTEGER(KIND=JPIM) :: ISTAT
CHARACTER(LEN=1) :: CLEC_PMON
ENERGY = 0
IF (MONINIT >= 0) THEN
   IF (MONINIT == 0) THEN ! The very first time only
      CALL GET_ENVIRONMENT_VARIABLE('EC_PMON',CLEC_PMON)
      IF (CLEC_PMON == '0') MONINIT = -2 ! Never try again
   ENDIF
   IF (MONINIT >= 0) THEN
      OPEN(503,FILE='/sys/cray/pm_counters/energy',IOSTAT=ISTAT,STATUS='old',ACTION='read')
      IF (ISTAT == 0) THEN
         READ(503,*,IOSTAT=ISTAT) ENERGY
         CLOSE(503)
         IF (ISTAT == 0) THEN
            IF (MONINIT == 0) THEN
               ENERGY_START = ENERGY
               MONINIT = 1 ! Ok
            ENDIF
            ENERGY = ENERGY - ENERGY_START
         ENDIF
      ENDIF
      IF (ISTAT /= 0) THEN
         MONINIT = -1 ! Never try again
         ENERGY = 0
      ENDIF
   ENDIF
ENDIF
POWER = 0
IF (MONINIT > 0) THEN
   OPEN(504,FILE='/sys/cray/pm_counters/power',IOSTAT=ISTAT,STATUS='old',ACTION='read')
   IF (ISTAT == 0) THEN
      READ(504,*,IOSTAT=ISTAT) POWER
      CLOSE(504)
   ENDIF
   IF (ISTAT /= 0) POWER = 0
ENDIF
END SUBROUTINE EC_PMON
