! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_MAP_SEND_PART1(YDIOS,KSIZEL,KSIZEG,YDFLDSC,KTAG,&
 & YDIOSMPP,KDOM_TYPE,PIOPROC1,PIOPROC2,PTSTEP)
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMIO_SERV_MAP_PLAN, ONLY : IO_SERV_SEND_PLAN
use ioflddesc_mod, only:&
 & ioflddesc
USE YOMIO_SERV, ONLY : IO_SERV
TYPE (IO_SERV), INTENT (INOUT) :: YDIOS
INTEGER (KIND=JPIM), INTENT (IN) :: KSIZEL
INTEGER (KIND=JPIM), INTENT (IN) :: KSIZEG
TYPE (IOFLDDESC), INTENT (IN) :: YDFLDSC (:)
INTEGER (KIND=JPIM), INTENT (IN) :: KTAG
TYPE (IO_SERV_SEND_PLAN), INTENT (INOUT) :: YDIOSMPP
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KDOM_TYPE
REAL (KIND=JPRB), INTENT (IN), OPTIONAL :: PIOPROC1
REAL (KIND=JPRB), INTENT (IN), OPTIONAL :: PIOPROC2
REAL (KIND=JPRB), INTENT (IN), OPTIONAL :: PTSTEP
call abor1("io_serv_map_send_part1.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_MAP_SEND_PART1
