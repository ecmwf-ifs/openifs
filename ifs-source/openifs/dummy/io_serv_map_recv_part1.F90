! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_MAP_RECV_PART1(YDIOS,KSIZEL,KSIZEG,YDFLDSC,PTSTEP,KTAG,YDIOSMPP,&
 & KDOM_TYPE,PIOPROC1,PIOPROC2,PTIME)
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMIO_SERV_MAP_PLAN, ONLY : IO_SERV_RECV_PLAN
use ioflddesc_mod , only:&
 & ioflddesc
USE YOMIO_SERV , ONLY : IO_SERV
TYPE (IO_SERV), INTENT (INOUT) :: YDIOS
INTEGER (KIND=JPIM), INTENT (IN) :: KSIZEL
INTEGER (KIND=JPIM), INTENT (IN) :: KSIZEG
TYPE (IOFLDDESC), INTENT (IN) :: YDFLDSC (:)
TYPE (IO_SERV_RECV_PLAN), INTENT (INOUT) :: YDIOSMPP
REAL (KIND=JPRB), INTENT (IN) :: PTSTEP
INTEGER (KIND=JPIM), INTENT (IN) :: KTAG
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KDOM_TYPE
REAL (KIND=JPRB), INTENT (IN), OPTIONAL :: PIOPROC1
REAL (KIND=JPRB), INTENT (IN), OPTIONAL :: PIOPROC2
REAL (KIND=JPRB), INTENT (IN), OPTIONAL :: PTIME
call abor1("io_serv_map_recv_part1.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_MAP_RECV_PART1
