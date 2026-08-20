! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUDIMO(YDGEOMETRY,KWINLEN,KULOUT,YDODB,CD_DBNAME)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jpim
USE DBASE_MOD, ONLY : DBASE
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KWINLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT
CLASS(DBASE) ,INTENT(OUT) :: YDODB
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: CD_DBNAME
call abor1("SUDIMO should never be called with this build configuration - EXITING")
END SUBROUTINE SUDIMO
