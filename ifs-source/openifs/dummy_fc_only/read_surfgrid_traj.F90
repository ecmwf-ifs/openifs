! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE READ_SURFGRID_TRAJ(YDGEOMETRY,YDSURF, YDDYN, YSURF, KSTEP,LDLASTRAJ,PTRAJ_BUF)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
use yomtraj , only:&
 & ngp5
USE YOMDYN , ONLY : TDYN
use surface_fields_mix , only:&
 & tsurf
USE ISO_C_BINDING
USE SRFSN_VGRID_MOD
USE SRFSN_REGRID_MOD
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KSTEP
LOGICAL ,INTENT(IN) :: LDLASTRAJ
TYPE(C_PTR) ,INTENT(IN) :: YSURF
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTRAJ_BUF(YDGEOMETRY%YRDIM%NPROMA,NGP5,YDGEOMETRY%YRDIM%NGPBLKS)
TYPE(TSURF) ,INTENT(INOUT) :: YDSURF
TYPE(TDYN) ,INTENT(INOUT) :: YDDYN
call abor1("READ_SURFGRID_TRAJ should never be called with this build configuration - EXITING")
END SUBROUTINE READ_SURFGRID_TRAJ
