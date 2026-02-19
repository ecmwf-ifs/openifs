! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CNT3_LAM(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDVARBC,KINITMONTH)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
USE MTRAJ_MOD , ONLY : MTRAJ
USE TYPE_MODEL , ONLY : MODEL
use parkind1 , only:&
 & jpim
USE VARBC_CLASS , ONLY : CLASS_VARBC
TYPE (GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
TYPE (FIELDS) ,INTENT(INOUT) :: YDFIELDS
TYPE (MTRAJ) ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (CLASS_VARBC), OPTIONAL, INTENT(INOUT) :: YDVARBC
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KINITMONTH
call abor1("cnt3_lam.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE CNT3_LAM
