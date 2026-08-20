! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE STORE_SEKF_CV(YDGEOMETRY,KIDIA,KFDIA,KLON,&
 & KSTGLO, KSTEP,&
 & PVSM_ML, PT2M_TEMP, PTD2M_TEMP,&
 & PTB_CMEM,&
 & KSTYPE)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 ,ONLY : JPIM, JPRB
use yomsekf , only:&
 & nslay
USE YOMSMOS , ONLY : NPOL_MAX, NANG_MAX
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB) ,INTENT(IN) :: PVSM_ML(YDGEOMETRY%YRDIM%NPROMA,NSLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: PT2M_TEMP(YDGEOMETRY%YRDIM%NPROMA), PTD2M_TEMP(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) ,INTENT(IN) :: PTB_CMEM(YDGEOMETRY%YRGEM%NGPTOT,0:NPOL_MAX-1,NANG_MAX)
INTEGER(KIND=JPIM),INTENT(IN) :: KSTYPE(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM),INTENT(IN) :: KLON, KIDIA, KFDIA, KSTEP
INTEGER(KIND=JPIM),INTENT(IN) :: KSTGLO
call abor1("STORE_SEKF_CV should never be called with this build configuration - EXITING")
END SUBROUTINE STORE_SEKF_CV
