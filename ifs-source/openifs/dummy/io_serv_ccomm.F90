! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IO_SERV_CCOMM (KRANK, KCOMMI, KCOMMO)
use parkind1, only:&
 & jpim
INTEGER(KIND=JPIM), INTENT (IN) :: KRANK (:)
INTEGER(KIND=JPIM), INTENT (IN) :: KCOMMI
INTEGER(KIND=JPIM), INTENT (OUT) :: KCOMMO
call abor1("io_serv_ccomm.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_CCOMM
