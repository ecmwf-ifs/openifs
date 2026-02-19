! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE EVDUVGEO(YDGEOMETRY,PVOR,PDIV,PMU,PMV,PSPU,PSPV)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jprb
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB) ,INTENT(IN) :: PVOR(:,:)
REAL(KIND=JPRB) ,INTENT(IN) :: PDIV(:,:)
REAL(KIND=JPRB) ,INTENT(IN) :: PMU(:)
REAL(KIND=JPRB) ,INTENT(IN) :: PMV(:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PSPU(:,:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PSPV(:,:)
call abor1("evduvgeo.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE EVDUVGEO
