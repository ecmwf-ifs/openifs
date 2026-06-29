! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!      -----------------------------------------------------------

! - Astronomical functions
! you will find the description in the annex 1 of the documentation
! RRS is the distance Sun-Earth
! RDS is the declination of the Earth
! RET is the equation of time

! Orbit of the earth

REAL(KIND=JPRD) :: RTETA,REL,REM,RRS,RLLS,RLLLS,RDS,RET
REAL(KIND=JPRD) :: PTIME,PTETA

RTETA(PTIME)=PTIME/(RDAY*365.25_JPRD)
REL(PTETA)=1.7535_JPRD+6.283076_JPRD*PTETA
REM(PTETA)=6.240075_JPRD+6.283020_JPRD*PTETA
RRS(PTETA)=REA*(1.0001_JPRD-0.0163_JPRD*SIN(REL(PTETA))&
          &+0.0037_JPRD*COS(REL(PTETA)))
! Relative movement Sun/Earth
RLLS(PTETA)=4.8951_JPRD+6.283076_JPRD*PTETA
RLLLS(PTETA)=4.8952_JPRD+6.283320_JPRD*PTETA-0.0075_JPRD*SIN(REL(PTETA))&
         &-0.0326_JPRD*COS(REL(PTETA))-0.0003_JPRD*SIN(2.0_JPRD*REL(PTETA))&
         &+0.0002_JPRD*COS(2.0_JPRD*REL(PTETA))
RDS(PTETA)=ASIN(SIN(REPSM)*SIN(RLLLS(PTETA)))
RET(PTETA)=591.8_JPRD*SIN(2.0_JPRD*RLLS(PTETA))-459.4_JPRD*SIN(REM(PTETA))&
  &+39.5_JPRD*SIN(REM(PTETA))*COS(2.0_JPRD*RLLS(PTETA))&
  &-12.7_JPRD*SIN(4._JPRD*RLLS(PTETA))-4.8_JPRD*SIN(2.0_JPRD*REM(PTETA))
!    -------------------------------------------------------------

