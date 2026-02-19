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

!      -----------------------------------------------------------

! - Astronomical functions
! you will find the description in the annex 1 of the documentation
! RRS is the distance Sun-Earth
! RDS is the declination of the Earth
! RET is the equation of time

! Orbit of the earth

REAL(KIND=JPRB) :: RTETA,REL,REM,RRS,RLLS,RLLLS,RDS,RET
REAL(KIND=JPRB) :: RRSAQUA, RDSAQUA, RETAQUA
REAL(KIND=JPRB) :: PTIME,PTETA

RTETA(PTIME)=PTIME/(RDAY*365.25_JPRB)
REL(PTETA)=1.7535_JPRB+6.283076_JPRB*PTETA
REM(PTETA)=6.240075_JPRB+6.283020_JPRB*PTETA
RRS(PTETA)=REA*(1.0001_JPRB-0.0163_JPRB*SIN(REL(PTETA))&
          &+0.0037_JPRB*COS(REL(PTETA)))

! Relative movement Sun/Earth
RLLS(PTETA)=4.8951_JPRB+6.283076_JPRB*PTETA
RLLLS(PTETA)=4.8952_JPRB+6.283320_JPRB*PTETA-0.0075_JPRB*SIN(REL(PTETA))&
         &-0.0326_JPRB*COS(REL(PTETA))-0.0003_JPRB*SIN(2.0_JPRB*REL(PTETA))&
         &+0.0002_JPRB*COS(2.0_JPRB*REL(PTETA))
RDS(PTETA)=ASIN(SIN(REPSM)*SIN(RLLLS(PTETA)))
RET(PTETA)=591.8_JPRB*SIN(2.0_JPRB*RLLS(PTETA))-459.4_JPRB*SIN(REM(PTETA))&
  &+39.5_JPRB*SIN(REM(PTETA))*COS(2.0_JPRB*RLLS(PTETA))&
  &-12.7_JPRB*SIN(4._JPRB*RLLS(PTETA))-4.8_JPRB*SIN(2.0_JPRB*REM(PTETA))

! APE aqua-planet special
RRSAQUA(PTETA)=REA
RDSAQUA(PTETA)=0.0_JPRB
RETAQUA(PTETA)=0.0_JPRB

!    -------------------------------------------------------------
