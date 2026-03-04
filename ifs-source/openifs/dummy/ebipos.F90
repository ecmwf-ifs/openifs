! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE EBIPOS(YDQTYPE,YDGEOMETRY,KFIELDS,LDBIP,PGPP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN
TYPE (TYPE_FPRQDYN), INTENT(IN) :: YDQTYPE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
LOGICAL, INTENT(IN) :: LDBIP(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PGPP(YDGEOMETRY%YRDIM%NPROMA,KFIELDS,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("ebipos.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE EBIPOS
