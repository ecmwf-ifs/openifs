! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE FP_SERV_INIT_PART1 (YDFPS, KCOMM_WORLD)
use parkind1, only:&
 & jpim
USE YOMFP_SERV, ONLY : FP_SERV
TYPE (FP_SERV), INTENT (INOUT) :: YDFPS
INTEGER (KIND=JPIM), INTENT (INOUT) :: KCOMM_WORLD
call abor1("fp_serv_init_part1.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE FP_SERV_INIT_PART1
