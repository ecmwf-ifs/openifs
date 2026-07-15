! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE READ_SURFGRID_TRAJ_FROMFA(YDGEOMETRY,YDSURF,YDMODEL,KSTEP,PTRAJ_SRFC)
USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
use surface_fields_mix , only:&
 & tsurf
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMTRAJ , ONLY : NGP5
TYPE(GEOMETRY) , INTENT(IN) :: YDGEOMETRY
TYPE(TSURF) , INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP
REAL(KIND=JPRB) , INTENT(OUT) :: PTRAJ_SRFC(YDGEOMETRY%YRDIM%NPROMA,NGP5,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("READ_SURFGRID_TRAJ_FROMFA should never be called with this build configuration - EXITING")
END SUBROUTINE READ_SURFGRID_TRAJ_FROMFA
