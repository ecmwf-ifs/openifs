! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

Subroutine IO_SERV_SUIOSCTMPL(YDIOS,YDGEOMETRY,PTSTEP,YDEWCOU,YDFPGEOMETRY,YDFPOPH)
USE GEOMETRY_MOD, ONLY : GEOMETRY
use parkind1, only:&
 & jprb
USE TYPE_FAOPH, ONLY : TFAOPH
USE YOMFPGEOMETRY, ONLY : TFPGEOMETRY
use yomio_serv, only:&
 & io_serv
USE GRIB_API_INTERFACE
USE YOEWCOU , ONLY : TEWCOU
TYPE (IO_SERV) , INTENT(INOUT) :: YDIOS
TYPE (GEOMETRY), INTENT(IN), OPTIONAL, TARGET :: YDGEOMETRY
REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PTSTEP
TYPE(TEWCOU), INTENT(INOUT), OPTIONAL :: YDEWCOU
TYPE (TFPGEOMETRY), INTENT(IN), OPTIONAL :: YDFPGEOMETRY
TYPE (TFAOPH), INTENT(IN), OPTIONAL :: YDFPOPH(:)
call abor1("io_serv_suiosctmpl.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE IO_SERV_SUIOSCTMPL
