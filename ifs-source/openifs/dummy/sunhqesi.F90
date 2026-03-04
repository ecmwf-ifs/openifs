! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUNHQESI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,KULOUT,LDEB)
USE YOMRIP , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jpim
USE YOMDYN , ONLY : TDYN
USE YOMDYNA , ONLY : TDYNA
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
TYPE(TDYN) ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA) ,INTENT(INOUT) :: YDDYNA
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT
LOGICAL ,INTENT(IN) :: LDEB
call abor1("sunhqesi.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUNHQESI
