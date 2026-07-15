! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUALLDFI(YDGEOMETRY,YDDIMF,KULOUT)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jpim
USE YOMDIMF , ONLY : TDIMF
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TDIMF) ,INTENT(INOUT) :: YDDIMF
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT
call abor1("SUALLDFI should never be called with this build configuration - EXITING")
END SUBROUTINE SUALLDFI
