! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE EUVGEOVD(YDGEOMETRY,PSPU,PSPV,PVOR,PDIV,PMU,PMV)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jprb
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB) ,INTENT(IN) :: PSPU(:,:)
REAL(KIND=JPRB) ,INTENT(IN) :: PSPV(:,:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PVOR(:,:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PDIV(:,:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PMU(:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PMV(:)
call abor1("euvgeovd.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE EUVGEOVD
