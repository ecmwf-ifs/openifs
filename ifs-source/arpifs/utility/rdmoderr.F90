! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE RDMODERR(YDMODERRCONF,YDGEOMETRY,YDRIP,CDFILE,YDSPCTLMODERR, &
    &               YDSPGPMODERR,YDSPERR,YDGPERR,CDFPATH, &
    &               LDRESTART,KNUM_PROCS_PER_IO_PROC)

!   Purpose.
!   --------
!     Read model error control variable in GRIB format.

!   Author.
!   -------
!     Y. Tremolet
!     Original   03-01-16

!   Modifications.
!   --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Tremolet    18-Mar-2004 Re-write
!        G. Radnoti    22-Nov-2006 Make YDGPERR optional
!        K. Yessad (July 2014): Move some variables.
!        M. Chrust (Jan 2020): OOPS cleaning
!        M. Chrust (Apr 2020): Add possibility to read restart files
! ------------------------------------------------------------------

USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE PARKIND1      , ONLY : JPRB, JPIM
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN        , ONLY : NULOUT
USE YOMRIP        , ONLY : TRIP
USE YOMMODERRCONF, ONLY : TMODERR_CONF
USE GRIDPOINT_FIELDS_MIX, ONLY : ASSIGNMENT(=), GRIDPOINT_FIELD
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD, SELF_MUL, SPECTRAL_NORM_LEVS, &
  &                             SPECTRAL_FIELD_CONTAINER

! ------------------------------------------------------------------

IMPLICIT NONE
TYPE(TMODERR_CONF),              INTENT(IN)    :: YDMODERRCONF
TYPE(GEOMETRY),                  INTENT(IN)    :: YDGEOMETRY
TYPE(TRIP),                      INTENT(INOUT) :: YDRIP
CHARACTER(LEN=*),                INTENT(IN)    :: CDFILE
TYPE(SPECTRAL_FIELD),  TARGET,   INTENT(INOUT) :: YDSPCTLMODERR
TYPE(SPECTRAL_FIELD),            INTENT(INOUT) :: YDSPGPMODERR
TYPE(SPECTRAL_FIELD),            INTENT(INOUT) :: YDSPERR
TYPE(GRIDPOINT_FIELD), OPTIONAL, INTENT(INOUT) :: YDGPERR
CHARACTER(LEN=*),      OPTIONAL, INTENT(IN)    :: CDFPATH
LOGICAL,               OPTIONAL, INTENT(IN)    :: LDRESTART
INTEGER(KIND=JPIM),    OPTIONAL, INTENT(IN)    :: KNUM_PROCS_PER_IO_PROC

! ------------------------------------------------------------------

LOGICAL :: LLRESTART
INTEGER(KIND=JPIM) :: NUM_PROCS_PER_IO_PROC
REAL(KIND=JPRB) :: ZH
TYPE(SPECTRAL_FIELD_CONTAINER) :: YLSPEC_CONTAINER(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ------------------------------------------------------------------

#include "read_spec_grib.intfb.h"
#include "read_spec_restart.intfb.h"
#include "spec_split.intfb.h"
#include "spec2grid.intfb.h"

! ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RDMODERR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
ASSOCIATE(TSTEP_ERR=>YDMODERRCONF%TSTEP_ERR,LSCALERR=>YDMODERRCONF%LSCALERR)
! ------------------------------------------------------------------

LLRESTART=.FALSE.
IF (PRESENT(LDRESTART)) LLRESTART=LDRESTART
NUM_PROCS_PER_IO_PROC=1
IF (PRESENT(KNUM_PROCS_PER_IO_PROC)) &
  & NUM_PROCS_PER_IO_PROC=KNUM_PROCS_PER_IO_PROC

! Read file
IF (LLRESTART) THEN
  YLSPEC_CONTAINER(1)%OBJ=>YDSPCTLMODERR
  IF(PRESENT(CDFPATH)) THEN
    CALL READ_SPEC_RESTART(YLSPEC_CONTAINER,CDFILE,CDFPATH=CDFPATH,&
      &                    KNUM_PROCS_PER_IO_PROC=NUM_PROCS_PER_IO_PROC)
  ELSE  
    CALL READ_SPEC_RESTART(YLSPEC_CONTAINER,CDFILE,&
      &                    KNUM_PROCS_PER_IO_PROC=NUM_PROCS_PER_IO_PROC)
  ENDIF  
ELSE
  CALL READ_SPEC_GRIB(YDGEOMETRY%YRMP,CDFILE,YDSPCTLMODERR)
ENDIF

! Model error is per hour in file. It might have to be converted
! to be per time-step.
IF (LSCALERR) THEN
  IF (TSTEP_ERR/=0.0_JPRB) THEN
    ZH = TSTEP_ERR/3600.0_JPRB
  ELSE
    ZH = TSTEP/3600.0_JPRB
  ENDIF
  CALL SELF_MUL(YDSPCTLMODERR,ZH)
ENDIF

! Print norms
CALL SPECTRAL_NORM_LEVS(YDSPCTLMODERR,'RDMODERR')

! Convert to proper SP/GP arrays
IF (PRESENT(YDGPERR)) THEN
  CALL SPEC_SPLIT(YDSPCTLMODERR,YDSPERR,YDSPGPMODERR)
  CALL SPEC2GRID(YDGEOMETRY,YDSPGPMODERR,YDGPERR)
ELSE
  YDSPERR=YDSPCTLMODERR
ENDIF

! ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RDMODERR',1,ZHOOK_HANDLE)
END SUBROUTINE RDMODERR
