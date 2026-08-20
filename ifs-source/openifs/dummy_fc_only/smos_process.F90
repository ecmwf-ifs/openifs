! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SMOS_PROCESS(YDGEM,PSMOS_OBS_BUF,KSMOS_BUF,PSMOS_ANGLE,&
 & PSMOS_FARAD,PSMOS_GEOMET)
 USE YOMGEM , ONLY : TGEM
 USE PARKIND1 ,ONLY : JPIM ,JPRB
 USE YOMDB
 use yomsmos , only:&
 & npol_max,&
 & nang_max
 TYPE(TGEM) ,INTENT(IN) :: YDGEM
 REAL(KIND=JPRB) ,INTENT(INOUT) :: PSMOS_OBS_BUF (YDGEM%NGPTOT, 0:NPOL_MAX-1, NANG_MAX)
 REAL(KIND=JPRB) ,INTENT(INOUT) :: PSMOS_ANGLE (YDGEM%NGPTOT, 0:NPOL_MAX-1, NANG_MAX)
 REAL(KIND=JPRB) ,INTENT(INOUT) :: PSMOS_FARAD (YDGEM%NGPTOT, 0:NPOL_MAX-1, NANG_MAX)
 REAL(KIND=JPRB) ,INTENT(INOUT) :: PSMOS_GEOMET (YDGEM%NGPTOT, 0:NPOL_MAX-1, NANG_MAX)
 INTEGER(KIND=JPIM) ,INTENT(INOUT) :: KSMOS_BUF (YDGEM%NGPTOT, 0:NPOL_MAX-1, NANG_MAX)
call abor1("SMOS_PROCESS should never be called with this build configuration - EXITING")
END SUBROUTINE SMOS_PROCESS
