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

! -- Field definitions for the radiation scheme driver --
!
! Life will be simpler when this can be an object that inherits from
! the base "field_definition" class
MODULE ECV_DEFINITIONS

USE PARKIND1              , ONLY: JPIM, JPRB
USE FIELD_DEFINITIONS_BASE, ONLY: SET_FVAR, TYPE_FVAR, FIELD_ACCESS_BASE, FIELD_METADATA_BASE
USE YOM_GRIB_CODES        , ONLY: NGRBALPHA, NGRBSKTECV, NGRBSSHECV, NGRBTSLECV,&
                                & NGRBSDFORE, NGRBFASGPPCOEF, NGRBFASRECCOEF

TYPE, EXTENDS(FIELD_ACCESS_BASE) :: ECVFIELD_ACCESS

  REAL(KIND=JPRB), POINTER :: ALPHA(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: SKTECV(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: SSHECV(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: TSLECV(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: SDFOR(:) => NULL()
  REAL(KIND=JPRB), POINTER :: GPPBFAS(:) => NULL()
  REAL(KIND=JPRB), POINTER :: RECBFAS(:) => NULL()

  CONTAINS

  PROCEDURE :: FIELD_MAP_STORAGE => ECVFIELD_MAP_STORAGE

END TYPE ECVFIELD_ACCESS

TYPE, EXTENDS(FIELD_METADATA_BASE) :: ECVFIELD_METADATA

  CONTAINS

  PROCEDURE :: FIELD_SET_METADATA => ECVFIELD_SET_METADATA
  PROCEDURE :: FIELD_GET_CLEVTYPE => ECV_FIELD_GET_CLEVTYPE

END TYPE ECVFIELD_METADATA

INTEGER(KIND=JPIM),PARAMETER :: JPNUMFIDS=41

TYPE TYPE_ECVFIELD_ID

  INTEGER(KIND=JPIM) :: ALPHA=10
  INTEGER(KIND=JPIM) :: SKTECV=20
  INTEGER(KIND=JPIM) :: SSHECV=21
  INTEGER(KIND=JPIM) :: TSLECV=22
  INTEGER(KIND=JPIM) :: SDFOR=30
  INTEGER(KIND=JPIM) :: GPPBFAS=31
  INTEGER(KIND=JPIM) :: RECBFAS=32

END TYPE TYPE_ECVFIELD_ID

TYPE(TYPE_ECVFIELD_ID), PARAMETER :: VID=TYPE_ECVFIELD_ID()

TYPE(ECVFIELD_METADATA) :: ECV_FIELD_METADATA

! 2nd dimension types
INTEGER(KIND=JPIM), PARAMETER :: VD2NONE         = 0  ! no second dimension
INTEGER(KIND=JPIM), PARAMETER :: VD2FULL         = 1 ! full model levels (third)
INTEGER(KIND=JPIM), PARAMETER :: VD2ALPHA_PERT   = 2 ! nb of alpha perturbetions (second dimension )
INTEGER(KIND=JPIM), PARAMETER :: VD2SKTSTEPS     = 3 ! Number of timesteps for skt
INTEGER(KIND=JPIM), PARAMETER :: VD2SSHSTEPS     = 4 ! Number of timesteps for ssh
INTEGER(KIND=JPIM), PARAMETER :: VD2TSLSTEPS     = 5 ! Number of timesteps for tsl

INTEGER(KIND=JPIM), PARAMETER :: NLEVELTYPES_ECV = 5

INTEGER(KIND=JPIM), PARAMETER :: NDIM3TYPES_ECV=1

#ifndef FIELD_MOD_TEST
#include "abor1.intfb.h"
#endif

!     ------------------------------------------------------------------

CONTAINS

!     ------------------------------------------------------------------

SUBROUTINE ECVFIELD_SET_METADATA(THIS, METADATA, KLEVELTYPES, KDIM3TYPES)
!
USE YOMHOOK          , ONLY : LHOOK, DR_HOOK, JPHOOK
!
CLASS(ECVFIELD_METADATA),     INTENT(IN)    :: THIS
TYPE(TYPE_FVAR), ALLOCATABLE, INTENT(INOUT) :: METADATA(:)
INTEGER(KIND=JPIM),           INTENT(  OUT) :: KLEVELTYPES
INTEGER(KIND=JPIM),           INTENT(  OUT) :: KDIM3TYPES
!
  call abor1("ECVFIELD_SET_METADATA in module ecv_definitions should never be called with OpenIFS - EXITING")

END SUBROUTINE ECVFIELD_SET_METADATA

!     ------------------------------------------------------------------

FUNCTION ECV_FIELD_GET_CLEVTYPE(THIS, KLEVELTYPE) RESULT(CLEVTYPE)
!
USE YOMHOOK          , ONLY : LHOOK, DR_HOOK, JPHOOK
!
CLASS(ECVFIELD_METADATA), INTENT(IN) :: THIS
INTEGER(KIND=JPIM)      , INTENT(IN) :: KLEVELTYPE
CHARACTER(LEN=3)  :: CLEVTYPE
  call abor1("ECV_FIELD_GET_CLEVTYPE (function) in module ecv_definitions should never be called with OpenIFS - EXITING")

END FUNCTION ECV_FIELD_GET_CLEVTYPE

! -------------------------------------------------
! Pointer mapping for the access type
! -------------------------------------------------

SUBROUTINE ECVFIELD_MAP_STORAGE(THIS, KID, STORAGE_1D, STORAGE_2D, STORAGE_3D, LD_NULLIFY)
!
USE YOMHOOK          , ONLY : LHOOK, DR_HOOK, JPHOOK
!
CLASS(ECVFIELD_ACCESS),            INTENT(INOUT) :: THIS
INTEGER(KIND=JPIM),                INTENT(IN)    :: KID
REAL(KIND=JPRB), OPTIONAL, TARGET, INTENT(IN)    :: STORAGE_1D(:), STORAGE_2D(:,:), STORAGE_3D(:,:,:)
LOGICAL, OPTIONAL,                 INTENT(IN)    :: LD_NULLIFY
  call abor1("ECVFIELD_MAP_STORAGE in module ecv_definitions should never be called with OpenIFS - EXITING")

END SUBROUTINE ECVFIELD_MAP_STORAGE

END MODULE ECV_DEFINITIONS
