! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE EFPFILTER(YDGEOMETRY,KFPDOM,LDFPFIL,KFMAX,PLTF,PFPFIL)
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE GEOMETRY_MOD , ONLY : GEOMETRY
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN) :: KFPDOM
LOGICAL, INTENT(IN) :: LDFPFIL(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KFMAX(KFPDOM)
REAL(KIND=JPRB), INTENT(IN) :: PLTF
REAL(KIND=JPRB), ALLOCATABLE, INTENT(OUT) :: PFPFIL(:,:)
call abor1("efpfilter.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE EFPFILTER
