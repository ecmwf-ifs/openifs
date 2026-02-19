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

MODULE ERRMOD_FIELDS_MOD
!> Model error fields (increment, state) for the IFS model
USE PARKIND1,             ONLY : JPIM, JPRB
USE YOMLUN,               ONLY : NULOUT, NULERR
USE TYPE_GEOMETRY,        ONLY : GEOMETRY
USE GEOMETRY_MOD,         ONLY : GEOMETRY_SAME
USE YOMMP0,               ONLY : MYPROC,NPROC
USE GRIDPOINT_FIELDS_MIX, ONLY : ASSIGNMENT(=), GRIDPOINT_FIELD
USE SPECTRAL_FIELDS_MOD,  ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE TYPE_MODEL,           ONLY : MODEL
USE YOMMODERRCONF,        ONLY : TMODERR_CONF
USE YOMHOOK,              ONLY : LHOOK, DR_HOOK, JPHOOK


IMPLICIT NONE
PRIVATE


TYPE,PUBLIC :: ERRMOD_FIELDS
  TYPE(TMODERR_CONF)                 :: CONF          ! Configuration
  TYPE(GRIDPOINT_FIELD), ALLOCATABLE :: GPMODERR(:)   ! Gridpoint model error
  TYPE(SPECTRAL_FIELD),  ALLOCATABLE :: SPMODERR(:)   ! Spectral model error
  TYPE(SPECTRAL_FIELD)               :: SPGPMODERR
  TYPE(SPECTRAL_FIELD)               :: SPCTLMODERR
  TYPE(GEOMETRY), POINTER            :: GEOM   => NULL()
  TYPE(MODEL), POINTER               :: MODEL  => NULL()
  LOGICAL                            :: LLAUX  = .FALSE.
  LOGICAL                            :: LLINCR = .TRUE.
END TYPE ERRMOD_FIELDS
END MODULE ERRMOD_FIELDS_MOD
