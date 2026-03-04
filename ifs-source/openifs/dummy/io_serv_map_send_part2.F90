! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_MAP_SEND_PART2 (YDIOS, YDIOSMPP)
USE YOMIO_SERV_MAP_PLAN, ONLY : IO_SERV_SEND_PLAN
USE YOMIO_SERV, ONLY : IO_SERV
TYPE (IO_SERV), INTENT (INOUT) :: YDIOS
TYPE (IO_SERV_SEND_PLAN), INTENT (INOUT) :: YDIOSMPP
call abor1("io_serv_map_send_part2.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_MAP_SEND_PART2
