! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_LOG (YDIOS, CDMESS, KTYPE, LDCLOSE)
use parkind1, only:&
 & jpim
use yomio_serv, only:&
 & io_serv
TYPE (IO_SERV), INTENT (INOUT) :: YDIOS
CHARACTER (LEN=*), INTENT (IN), OPTIONAL :: CDMESS
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KTYPE
LOGICAL, INTENT (IN), OPTIONAL :: LDCLOSE
call abor1("io_serv_log.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_LOG
