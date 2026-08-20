! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SMAP_PROCESS(YDGEM,PSMAP_OBS_BUF,KSMAP_BUF,PSMAP_ANGLE)
 USE YOMGEM , ONLY : TGEM
 USE PARKIND1 ,ONLY : JPIM, JPRB
 USE YOMDB
 use yomsmos , only:&
 & npol_max
 TYPE(TGEM) ,INTENT(IN) :: YDGEM
 REAL(KIND=JPRB) ,INTENT(INOUT) :: PSMAP_OBS_BUF (YDGEM%NGPTOT, 0:NPOL_MAX-1)
 REAL(KIND=JPRB) ,INTENT(INOUT) :: PSMAP_ANGLE (YDGEM%NGPTOT, 0:NPOL_MAX-1)
 INTEGER(KIND=JPIM) ,INTENT(INOUT) :: KSMAP_BUF (YDGEM%NGPTOT, 0:NPOL_MAX-1)
call abor1("SMAP_PROCESS should never be called with this build configuration - EXITING")
END SUBROUTINE SMAP_PROCESS
