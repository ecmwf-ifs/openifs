! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_RECV_REQ (YDIOS, YDREQ, KSTEPMAX)
use parkind1, only:&
 & jpim
USE YOMIO_SERV, ONLY : IO_SERV
USE YOMIO_SERV_REQ, ONLY : IO_SERV_REQ
TYPE (IO_SERV), INTENT (INOUT) :: YDIOS
TYPE (IO_SERV_REQ), POINTER :: YDREQ
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: KSTEPMAX
call abor1("io_serv_recv_req.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_RECV_REQ
