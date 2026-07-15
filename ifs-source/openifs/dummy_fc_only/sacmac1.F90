! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SACMAC1(YDDIMF,YDODB,YDGEOMETRY,PFJPCOST,YDJOT)
USE YOMDIMF , ONLY : TDIMF
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jprb
USE JO_TABLE_MOD , ONLY : JO_TABLE
USE DBASE_MOD, ONLY : DBASE
TYPE(TDIMF) ,INTENT(INOUT) :: YDDIMF
CLASS(DBASE) ,INTENT(INOUT) :: YDODB
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB) ,INTENT(IN) :: PFJPCOST
TYPE(JO_TABLE) ,INTENT(INOUT) :: YDJOT
call abor1("SACMAC1 should never be called with this build configuration - EXITING")
END SUBROUTINE SACMAC1
