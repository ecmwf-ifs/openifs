! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_INIT_PART2 (YDIOS, KCOMM_WORLD)
USE YOMIO_SERV, ONLY : IO_SERV
use parkind1, only:&
 & jpim
TYPE (IO_SERV), INTENT (INOUT) :: YDIOS
INTEGER (KIND=JPIM), INTENT (INOUT) :: KCOMM_WORLD
call abor1("io_serv_init_part2.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_INIT_PART2
