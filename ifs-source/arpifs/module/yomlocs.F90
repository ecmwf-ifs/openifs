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

MODULE YOMLOCS

! Manage observation location information that is extracted from the observation
! object (from the ODB) and required to create the model-to-obs interpolation
! object (the GOMs).

USE PARKIND1    , ONLY : JPIM, JPRB, JPRD
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN      , ONLY : NULOUT, NULERR, NULNAM
USE YOMCT0      , ONLY : LALLOPR, LCANARI, LELAM, LRPLANE
USE YOMANCS     , ONLY : RMDI
USE YOMCST      , ONLY : RPI
USE YOMCOCTP    , ONLY : NSATEM, NALLSKY, NLIMB
USE YOMSATS     , ONLY : SATGRP_TABLE
USE GOM2D_SUPPORT, ONLY: JP_TREAT_2DGOM_SINGLE, JP_TREAT_2DGOM_LIMB_PLANE, JP_TREAT_2DGOM_SLANT_PATH,&
 &                       JP_TREAT_2DGOM_UNSPEC, LIMB_PLANE
USE YOMMP0      , ONLY : NPRINTLEV, MYPROC, NPROC
! for radar reflectivity
USE YOMECTAB    , ONLY : NTSLPT, NTSLTOBP, NOBSTYP
USE PARDIMO     , ONLY : JPNOTP
USE IFS_DBASE_VIEW_MOD, ONLY : IFS_DBASE_VIEW
USE DBASE_MOD   , ONLY : DBASE

IMPLICIT NONE





SAVE

TYPE TLOCS
  REAL(KIND=JPRB)    :: LAT      ! lat@hdr
  REAL(KIND=JPRB)    :: LON      ! lon@hdr
  INTEGER(KIND=JPIM) :: SEQNO    ! seqno@hdr
  INTEGER(KIND=JPIM) :: OBSTYPE  ! obstype@hdr
  INTEGER(KIND=JPIM) :: CODETYPE ! codetype@hdr
  INTEGER(KIND=JPIM) :: AREATYPE ! areatype@hdr
  INTEGER(KIND=JPIM) :: RETRTYPE ! retrtype@hdr
  INTEGER(KIND=JPIM) :: OBSPROF  ! from MKGLOBSTAB_OB
  INTEGER(KIND=JPIM) :: PROFNUM  ! from MKGLOBSTAB_OB
END TYPE TLOCS

!------------------------------------------------------------------------
! Full set of observation location info in form suitable for input to
! the model-obs interpolation operator (GOMs)
TYPE :: TOBSLOCS

  INTEGER(KIND=JPIM), POINTER:: ISTAPROFS(:)=>NULL(), ITSLTPROF(:)=>NULL()
  TYPE(TLOCS), POINTER :: YLOCS(:)=>NULL()
  INTEGER(KIND=JPIM), POINTER :: IOTP(:)=>NULL(),INPROFS(:)=>NULL()
  INTEGER(KIND=JPIM) :: NACTIM            ! number of timeslots
  INTEGER(KIND=JPIM) :: NOBPROF           ! NUMBER OF OBS.EQ. PROFILES ON THIS PE
  INTEGER(KIND=JPIM) :: NOBPROFG          ! GLOBAL NUMBER OF OBS.EQ. PROFILES
  INTEGER(KIND=JPIM) :: NOBSPROFS(JPNOTP) ! MAXIMUM NUMBER OF OBS.EQ. PROFILES PER OBS. TYPE

END TYPE
  ! FORECAST_ONLY or OIFS - removed obs routines to enable build in FC / OIFS mode


END MODULE YOMLOCS
