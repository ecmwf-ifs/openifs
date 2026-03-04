! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE FP_SERV_CPFPFILTER (YDFPS, YDGEOMETRY, YDFPFILTERS)
USE YOMFP_SERV, ONLY : FP_SERV
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMFPFILTERS, ONLY : TFPFILTERS
TYPE (FP_SERV), INTENT (INOUT), TARGET :: YDFPS
TYPE (GEOMETRY), INTENT (IN) :: YDGEOMETRY
TYPE (TFPFILTERS),INTENT (INOUT) :: YDFPFILTERS
call abor1("fp_serv_cpfpfilter.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE FP_SERV_CPFPFILTER
