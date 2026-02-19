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

! Prototype for managing model fields in the IFS
!
! Much has been copied from (or inspired by) the new GOMs  - obviously
! the GOMs would eventually use the field module stuff, not duplicate it!
!
! This code is far from perfect yet and there is some scope for performance
! tuning, but that can wait until profiling tells us what is actually slow.
!
! AJG 24/10/2012
!
MODULE FIELD_CONTAINER_GP_MOD

USE PARKIND1              , ONLY : JPIM, JPRB
USE YOMHOOK               , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN                , ONLY : NULOUT
#ifndef FIELD_CONTAINER_GP_MOD_TEST
USE YOMGFL                , ONLY : TGFL
USE YOMGMV                , ONLY : TGMV
USE SURFACE_FIELDS_MIX    , ONLY : TSURF
#endif
USE FIELD_DEFINITIONS_BASE, ONLY : TYPE_FVAR,TYPE_CONTROL,JP_T0_DL,JP_T0_DM,JP_T9_DL,JP_T9_DM,JP_T0,JP_T9,&
 &                                 FIELD_ACCESS_BASE, FIELD_METADATA_BASE
USE FIELD_DEFINITIONS     , ONLY : FIELD_ACCESS
#ifndef FIELD_CONTAINER_GP_MOD_TEST
USE FIELD_GFL_WRAPPER     , ONLY : GFL_MAPPING_SETUP, GFL_FID_SETUP, GFL_ATTACH, GFL_AS_A_WHOLE,GFL_MAP
#endif
USE FIELD_CONTAINER_BASE_MOD

IMPLICIT NONE

PRIVATE

! -------------------------------------------------
! Public interface
! -------------------------------------------------


! -------------------------------------------------


INTEGER(KIND=JPIM), PARAMETER :: JP_ALL_BLOCKS = 0

! Storage scheme: NFTBC eventually, delete all but the best one.
!INTEGER(KIND=JPIM), PARAMETER :: JP_STORAGE1D  = 1 ! basic 1d storage scheme (not rank-conformant at f95)
INTEGER(KIND=JPIM), PARAMETER :: JP_STORAGE3D  = 4 ! "3d" storage scheme preserving rank-conformance
!INTEGER(KIND=JPIM), PARAMETER :: JP_STORAGEOLD = 3 ! backward compatibility: use old gfl, gmvs etc.

INTEGER(KIND=JPIM), PARAMETER :: JP_FADD=1,JP_FCOPY=2
! 5) Derived types, constants
PUBLIC :: FIELD_CONTAINER_GP, FIELD_CONTAINER_GP_CONTAINER

TYPE,EXTENDS(TYPE_FIELD_INDEX_BASE) ::  TYPE_FIELD_INDEX

  ! For one variable
  INTEGER(KIND=JPIM) :: ISTART=-1
  INTEGER(KIND=JPIM) :: IEND=-1
  LOGICAL :: L_HDERS = .FALSE.
  TYPE(TYPE_FIELD_INDEX),POINTER :: DM !=> NULL()
  TYPE(TYPE_FIELD_INDEX),POINTER :: DL => NULL()

END TYPE TYPE_FIELD_INDEX

TYPE,EXTENDS(FIELD_CONTAINER_BASE) :: FIELD_CONTAINER_GP
  PRIVATE

  ! Storage type (e.g. 1D, 3D, attach to old GFLS etc.)
#ifndef FIELD_CONTAINER_GP_MOD_TEST
  TYPE(TGMV),POINTER         :: YRGMV=>NULL()
  TYPE(TGFL),POINTER         :: YRGFL=>NULL()
  TYPE(TSURF),POINTER        :: YRSURF=>NULL()
#endif

  ! 1D Storage area and its size; indexes into the storage area
  INTEGER(KIND=JPIM)             :: ILEN_STORAGE1D
  REAL(KIND=JPRB), ALLOCATABLE   :: STORAGE1D(:)
  TYPE(TYPE_FIELD_INDEX),ALLOCATABLE :: INDEX_GP(:)
  TYPE(TYPE_FIELD_INDEX),ALLOCATABLE :: INDEX_BLOCK(:,:)

  ! GFL-mapping storage area NFTBC
  LOGICAL :: L_ATTACH_GFL_NPROMA_FIXED  ! TRUE = turns OFF the last block hack
  REAL(KIND=JPRB), ALLOCATABLE :: STORAGE_LAST_BLOCK(:,:,:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: INDEX_LAST_BLOCK(:)
  REAL(KIND=JPRB), ALLOCATABLE    :: UNASSIGNED_GFL(:,:,:) ! to deal with unused fields

  ! Properties of all fields that have been defined
  INTEGER(KIND=JPIM)          :: NHDERS=0

  ! Basic dimensions
  INTEGER(KIND=JPIM)          :: NBLOCKS
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPOINTS_BLOCK(:) ! "nproma" kind of thing

  ! List of all the fields that actually in use
  LOGICAL                     :: LHDERS

  CONTAINS
! 1) Creator and destructor
  PROCEDURE ::  FIELD_CREATE
  PROCEDURE ::  FIELD_DESTROY
! 2) Direct field access
  PROCEDURE :: FIELD_OPEN_FAC
  PROCEDURE :: FIELD_CLOSE_FAC
  PROCEDURE :: FIELD_OPEN_ONE_1D
  PROCEDURE :: FIELD_OPEN_ONE_2D
  PROCEDURE :: FIELD_OPEN_ONE_3D
  GENERIC   :: FIELD_OPEN_ONE => FIELD_OPEN_ONE_1D
  GENERIC   :: FIELD_OPEN_ONE => FIELD_OPEN_ONE_2D
  GENERIC   :: FIELD_OPEN_ONE => FIELD_OPEN_ONE_3D
  PROCEDURE :: FIELD_CLOSE_ONE_1D
  PROCEDURE :: FIELD_CLOSE_ONE_2D
  PROCEDURE :: FIELD_CLOSE_ONE_3D
  GENERIC   :: FIELD_CLOSE_ONE => FIELD_CLOSE_ONE_1D
  GENERIC   :: FIELD_CLOSE_ONE => FIELD_CLOSE_ONE_2D
  GENERIC   :: FIELD_CLOSE_ONE => FIELD_CLOSE_ONE_3D
  PROCEDURE :: FIELD_OPEN_MULTI
  PROCEDURE :: FIELD_CLOSE_MULTI
  PROCEDURE :: FIELD_ITERATED
! 3) Field operators (indirect access)
  PROCEDURE :: FIELD_EQUALS
  PROCEDURE :: FIELD_ADD
  PROCEDURE :: FIELD_COPY
  PROCEDURE :: FIELD_NORM
! 4) Multi-threading helpers; inquiry functions
  PROCEDURE :: FIELD_BLOCK_SIZE
  PROCEDURE :: FIELD_NBLOCKS
!  PROCEDURE :: FIELD_ACTIVE_CONTAINER
  PROCEDURE :: KINDEX_OLD
END TYPE FIELD_CONTAINER_GP

TYPE FIELD_CONTAINER_GP_CONTAINER
  CLASS(FIELD_CONTAINER_GP), POINTER :: OBJ
END TYPE FIELD_CONTAINER_GP_CONTAINER    


#ifndef FIELD_CONTAINER_GP_MOD_TEST
#include "abor1.intfb.h"
#endif

CONTAINS

SUBROUTINE FIELD_CREATE(SELF,NAMESPACE,KPOINTS,KLEVELS,KPROMA,KFIDS,KFSET,KATTACH,YDGMV,YDGFL,YDSURF,&
  & LD_ATTACH_GFL_NPROMA_FIXED,LDHDERS,KDIM3S)

! -------------------------------------------------
! Create a field container
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP) ,   INTENT(INOUT) :: SELF
CLASS(FIELD_METADATA_BASE),   INTENT(IN)    :: NAMESPACE  ! Metadata defining what namespace we're using (Main, radiation, etc.)
INTEGER(KIND=JPIM),           INTENT(IN)    :: KPOINTS    ! horizontal size of field on this PE
INTEGER(KIND=JPIM),           INTENT(IN)    :: KLEVELS(:) ! vertical dimensions
INTEGER(KIND=JPIM),           INTENT(IN)    :: KPROMA     ! size of blocks for multi-threading
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KFIDS(:)    ! list of field IDs (otherwise make all fields)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KFSET       ! field definition set
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KATTACH     ! Attachment mode (GFL backwards compatibility)
#ifndef FIELD_CONTAINER_GP_MOD_TEST
TYPE(TGMV),OPTIONAL,TARGET ,INTENT(IN)      :: YDGMV
TYPE(TGFL) ,OPTIONAL,TARGET  ,INTENT(IN)    :: YDGFL
TYPE(TSURF) ,OPTIONAL,TARGET  ,INTENT(IN)   :: YDSURF
#else
LOGICAL, OPTIONAL, INTENT(INOUT) :: YDGMV  ! Just to avoid having to mess with dummy argument list above
LOGICAL, OPTIONAL, INTENT(INOUT) :: YDGFL  ! Just to avoid having to mess with dummy argument list above
LOGICAL, OPTIONAL, INTENT(INOUT) :: YDGMV  ! Just to avoid having to mess with dummy argument list above
#endif
LOGICAL           , OPTIONAL, INTENT(IN)    :: LD_ATTACH_GFL_NPROMA_FIXED ! Do not adjust last nproma size
LOGICAL           , OPTIONAL, INTENT(IN)    :: LDHDERS ! Include horizontal derivatives
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KDIM3S(:)

INTEGER(KIND=JPIM) :: ID, ISIZE_ALL, ISIZE_VAR, IBLOCK, INLEVELS, I, I2D
INTEGER(KIND=JPIM) :: ISIZE_BLOCK, IMAXLEV
INTEGER(KIND=JPIM) :: IKGLO
LOGICAL,ALLOCATABLE :: LLON(:),LLACTIVE(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CREATE',0,ZHOOK_HANDLE)

! Define namespace and metadata
! Basic dimensions
CALL SELF%FIELD_CREATE_BASE(NAMESPACE,KPOINTS,KLEVELS,KFIDS,KATTACH,KSTORAGE_TYPE=JP_STORAGE_ARR,KDIM3S=KDIM3S)

SELF%NHDERS = 0
SELF%LHDERS = .FALSE.
IF(PRESENT(LDHDERS)) THEN
  SELF%LHDERS = LDHDERS
ENDIF

! Set up the "nproma" parallelisation
SELF%NBLOCKS = (SELF%NPOINTS-1)/KPROMA+1
ALLOCATE(SELF%NPOINTS_BLOCK(SELF%NBLOCKS))
SELF%NPOINTS_BLOCK(:) = KPROMA
IF (MOD(SELF%NPOINTS,SELF%NBLOCKS) /= 0) THEN
  IKGLO = 1+(SELF%NBLOCKS-1)*KPROMA
  SELF%NPOINTS_BLOCK(SELF%NBLOCKS) = MIN(KPROMA,SELF%NPOINTS-IKGLO+1)
ENDIF

! Initialise the index

ALLOCATE(SELF%INDEX_GP(SELF%NFIELDS))


! Initiliase the locking and iterator mechanism

! Storage type and its initialisation
IF (PRESENT(KATTACH)) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

  ! For initial implementation in the IFS, we will
  ! point at the old GFL/GMV/GMVS anyway, to allow the progressive
  ! implementation of the new fields in the IFS.
  SELF%YRGMV => YDGMV
  SELF%YRGFL => YDGFL
  SELF%YRSURF => YDSURF
  ! Set up the mappings onto the old GFL (etc.) pointers
  CALL GFL_MAPPING_SETUP(SELF%GFL_MAP_INTERNAL,SELF%YRGMV,SELF%YRGFL,SELF%YRGFL%YGFL,SELF%YRSURF,SELF%METADATA)

  ! Allow the existing GFL/GMV/surface arrays to define what's on or off
  ALLOCATE(LLON(SELF%NFIELDS))
  ALLOCATE(LLACTIVE(SELF%NFIELDS))
  DO ID =1,SELF%NFIELDS
    LLON(ID) = SELF%FIELD_INDEX(ID)%L_ON
    LLACTIVE(ID) = SELF%FIELD_INDEX(ID)%L_ACTIVE
  ENDDO

  CALL GFL_FID_SETUP(SELF%GFL_MAP_INTERNAL,SELF%NFIELDS, KATTACH, SELF%NFIELDS_ON, LLON, LLACTIVE,SELF%MALL_FIELDS, &
   & KFIDS)
  DO ID =1,SELF%NFIELDS
    SELF%FIELD_INDEX(ID)%L_ON = LLON(ID)
    SELF%FIELD_INDEX(ID)%L_ACTIVE = LLACTIVE(ID)
    IF(SELF%LHDERS)  THEN
      SELF%INDEX_GP(ID)%L_HDERS = SELF%FIELD_INDEX(ID)%L_ON
    ELSE
      SELF%INDEX_GP(ID)%L_HDERS = .FALSE.
    ENDIF
  ENDDO

  SELF%M_ATTACH   = KATTACH
  SELF%L_ATTACH_GFL_NPROMA_FIXED = .FALSE.
  IF(PRESENT(LD_ATTACH_GFL_NPROMA_FIXED)) &
    & SELF%L_ATTACH_GFL_NPROMA_FIXED = LD_ATTACH_GFL_NPROMA_FIXED

#else

  WRITE(0,*) 'ATTACHING TO GFLS IS NOT SUPPORTED'
  CALL ABORT

#endif

ELSE

  ! Ultimately, we will probably use a 1D storage area with F2003 rank-compliant
  ! pointer mappings. There is also a possibility to use a series of 3D storage
  ! areas in the future, which will preserve rank-compliance
  SELF%L_ATTACH_GFL_NPROMA_FIXED = .FALSE.

  IF(PRESENT(KFIDS)) THEN
    ! Turn on the required fields.
    ! NFTBC Here, fids DEFINES what is on and what is off (unlike GFL/GMV approach)
    SELF%NFIELDS_ON = SIZE(KFIDS)
    DO I=1,SELF%NFIELDS_ON
      ID = KFIDS(I)
      SELF%FIELD_INDEX(ID)%L_ON = .TRUE.
      SELF%FIELD_INDEX(ID)%L_ACTIVE = .TRUE.
      IF(SELF%LHDERS)  SELF%INDEX_GP(ID)%L_HDERS = .TRUE.
    ENDDO
  ELSE
    ! Turn on all fields
    SELF%NFIELDS_ON = SELF%NFIELDS
    SELF%FIELD_INDEX(:)%L_ON = .TRUE.
    SELF%FIELD_INDEX(:)%L_ACTIVE = .TRUE.
    IF(SELF%LHDERS)  SELF%INDEX_GP(:)%L_HDERS = .TRUE.
  ENDIF

  ! Keep a list of the fields that are turned on
  ALLOCATE(SELF%MALL_FIELDS(SELF%NFIELDS_ON))
  IF(PRESENT(KFIDS)) THEN
    SELF%MALL_FIELDS = KFIDS
  ELSE
    DO ID=1,SELF%NFIELDS
      SELF%MALL_FIELDS(ID) = ID
    ENDDO
  ENDIF

ENDIF

IF(SELF%MSTORTYPE == JP_STORAGE1D) THEN

  ! Set up the index into the 1D storage area
  !
  ! The 1D approach means that fields can be any shape, e.g. 1D, 2D, and
  ! there is no wasted space when npoints is not an exact multiple of nproma.

  ! Initialise the by-block index
  ALLOCATE(SELF%INDEX_BLOCK(SELF%NFIELDS,SELF%NBLOCKS))
  ISIZE_ALL = 0
  DO IBLOCK = 1,SELF%NBLOCKS
    DO ID=1,SELF%NFIELDS
      IF (SELF%FIELD_INDEX(ID)%L_ON) THEN

        IF (SELF%METADATA(ID)%NDIMS == 1) THEN
          INLEVELS = 1
        ELSEIF (SELF%METADATA(ID)%NDIMS == 2) THEN
          INLEVELS = SELF%NLEVELS(SELF%METADATA(ID)%D2TYPE)
        ELSEIF (SELF%METADATA(ID)%NDIMS == 3) THEN
          INLEVELS = SELF%NLEVELS(SELF%METADATA(ID)%D2TYPE) * SELF%NDIM3S(SELF%METADATA(ID)%D3TYPE)
        ENDIF

        ISIZE_VAR = SELF%NPOINTS_BLOCK(IBLOCK) * INLEVELS

        SELF%INDEX_BLOCK(ID,IBLOCK)%ISTART = ISIZE_ALL + 1
        SELF%INDEX_BLOCK(ID,IBLOCK)%IEND   = ISIZE_ALL + ISIZE_VAR
        ISIZE_ALL                       = ISIZE_ALL + ISIZE_VAR
        IF(SELF%INDEX_GP(ID)%L_HDERS) THEN
          ALLOCATE(SELF%INDEX_BLOCK(ID,IBLOCK)%DL)
          SELF%INDEX_BLOCK(ID,IBLOCK)%DL%ISTART = ISIZE_ALL + 1
          SELF%INDEX_BLOCK(ID,IBLOCK)%DL%IEND   = ISIZE_ALL + ISIZE_VAR
          ISIZE_ALL                          = ISIZE_ALL + ISIZE_VAR
          ALLOCATE(SELF%INDEX_BLOCK(ID,IBLOCK)%DM)
          SELF%INDEX_BLOCK(ID,IBLOCK)%DM%ISTART = ISIZE_ALL + 1
          SELF%INDEX_BLOCK(ID,IBLOCK)%DM%IEND   = ISIZE_ALL + ISIZE_VAR
          ISIZE_ALL                          = ISIZE_ALL + ISIZE_VAR
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  SELF%ILEN_STORAGE1D = ISIZE_ALL
  ALLOCATE(SELF%STORAGE1D(SELF%ILEN_STORAGE1D))
ELSEIF(SELF%MSTORTYPE == JP_STORAGE_ARR) THEN
  CALL SELF%ALLOCATE_STORAGE_ARRAYS(KPROMA=KPROMA)
ELSEIF(SELF%MSTORTYPE == JP_STORAGEOLD) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

  ! Transitional arrangements: just act as a wrapper to GFL/GMV/surface

  ! Most pointers can be directly mapped into the GFL, but to deal
  ! with the "last" nproma block and the non-rank-compliant
  ! nature of much of the IFS, we need to create non-strided
  ! array storage for the last nproma block of all 2D variables.
  ! ASSUMPTIONS: - all 2D GFL-style arrays have ==NFLEVG levels
  !              - last NPROMA block size is globally constant
  I2D         = COUNT(SELF%FIELD_INDEX(:)%L_ON .AND. SELF%METADATA(:)%NDIMS == 2)
  ISIZE_BLOCK = MINVAL(SELF%NPOINTS_BLOCK)
  IMAXLEV     = SELF%NLEVELS(1) ! AJGDB/NFTBC: Hidden field_definitions dependence - must be NFLEVG!
  ALLOCATE(SELF%STORAGE_LAST_BLOCK(ISIZE_BLOCK,IMAXLEV,I2D))
  ALLOCATE(SELF%INDEX_LAST_BLOCK(SELF%NFIELDS))

  ! Index that will translate a FID into a position in the last block storage
  I2D = 1
  DO ID=1,SELF%NFIELDS
    IF(SELF%FIELD_INDEX(ID)%L_ON .AND. SELF%METADATA(ID)%NDIMS == 2) THEN
      SELF%INDEX_LAST_BLOCK(ID) = I2D
      I2D = I2D+1
    ELSE
      SELF%INDEX_LAST_BLOCK(ID) = -1
    ENDIF
  ENDDO

  ! A solution to the "unused GFL" problem. For unassigned GFLs
  ! this allows us to pass pointers to unused GFL fields through
  ! F77 interfaces, where a nullified pointer would cause a crash.
  ALLOCATE(SELF%UNASSIGNED_GFL( MAXVAL(SELF%NPOINTS_BLOCK),IMAXLEV,1))
  SELF%UNASSIGNED_GFL = -HUGE(SELF%UNASSIGNED_GFL)

#endif

ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CREATE',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CREATE
!=======================================================================================
SUBROUTINE FIELD_DESTROY(SELF)
! -------------------------------------------------
! Destroy a field container
! -------------------------------------------------
CLASS(FIELD_CONTAINER_GP), INTENT(INOUT) :: SELF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_DESTROY',0,ZHOOK_HANDLE)

IF(SELF%LCREATED) THEN
  IF(SELF%MSTORTYPE == JP_STORAGE1D) THEN
    DEALLOCATE(SELF%STORAGE1D)
    DEALLOCATE(SELF%INDEX_BLOCK)
  ELSEIF(SELF%MSTORTYPE == JP_STORAGE_ARR) THEN
    CALL SELF%DEALLOCATE_STORAGE_ARRAYS()
  ELSEIF(SELF%MSTORTYPE == JP_STORAGEOLD) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

    DEALLOCATE(SELF%STORAGE_LAST_BLOCK, SELF%INDEX_LAST_BLOCK)

  ! Check to see that unassigned GFLs have not been written to.
    IF(ANY(SELF%UNASSIGNED_GFL /= -HUGE(SELF%UNASSIGNED_GFL))) &
     & CALL ABOR1('FIELD_CONTAINER_GP_MOD: TRANSITIONAL GFL MAPPING TO UNASSIGNED GFL DOES NOT PERMIT WRITE ACCESS')

    DEALLOCATE(SELF%UNASSIGNED_GFL)

#endif

  ENDIF
  DEALLOCATE(SELF%INDEX_GP,SELF%MALL_FIELDS,SELF%NLEVELS,SELF%METADATA,&
   & SELF%NPOINTS_BLOCK,SELF%FIELD_INDEX)
  SELF%LCREATED = .FALSE.
ENDIF

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_DESTROY',1,ZHOOK_HANDLE)
END SUBROUTINE FIELD_DESTROY
!=======================================================================================
SUBROUTINE FIELD_OPEN_FAC(SELF, KFIDS, KBLOCK, FAC, KBLOCK_SIZE,KFIDS_HD)
! -------------------------------------------------
! Get access to model fields via an access type
!   - There is one wrapper procedure for each access type
!   - The common parts are found in field_open_common
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),           INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM), TARGET, OPTIONAL,INTENT(IN)    :: KFIDS(:)
INTEGER(KIND=JPIM),                  INTENT(IN)    :: KBLOCK
CLASS(FIELD_ACCESS_BASE),            INTENT(OUT)   :: FAC
INTEGER(KIND=JPIM), OPTIONAL,        INTENT(OUT)   :: KBLOCK_SIZE
INTEGER(KIND=JPIM), TARGET, OPTIONAL,INTENT(IN)    :: KFIDS_HD(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_FAC',0,ZHOOK_HANDLE)

SELECT TYPE (FAC)
TYPE IS (FIELD_ACCESS)
  IF(SELF%LHDERS.AND.PRESENT(KFIDS_HD)) THEN
    IF(.NOT. ASSOCIATED(FAC%DL)) ALLOCATE(FIELD_ACCESS::FAC%DL)
    IF(.NOT. ASSOCIATED(FAC%DM)) ALLOCATE(FIELD_ACCESS::FAC%DM)
  ENDIF
END SELECT

CALL FIELD_OPEN_COMMON(SELF, FAC%ACCESS_CONTROL, FAC%FIELDS, KFIDS=KFIDS, KBLOCK=KBLOCK, KBLOCK_SIZE=KBLOCK_SIZE)
CALL MAP_NAMES_TO_STORAGE(SELF, FAC%FIELDS, KBLOCK, MFAC=FAC,KFIDS_HD=KFIDS_HD)

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_FAC',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_FAC
!=======================================================================================
SUBROUTINE FIELD_OPEN_COMMON(SELF, ACCESS_CONTROL, KFIELDS, KFIDS, KBLOCK, KBLOCK_SIZE)

CLASS(FIELD_CONTAINER_GP),            INTENT(INOUT) :: SELF
TYPE(TYPE_CONTROL),                   INTENT(OUT)   :: ACCESS_CONTROL
INTEGER(KIND=JPIM), POINTER          ,INTENT(OUT)   :: KFIELDS(:) ! OUT: a list of fields that will be available in the access type
INTEGER(KIND=JPIM), OPTIONAL, TARGET, INTENT(IN)    :: KFIDS(:) ! list of field IDs (otherwise get all fields)
INTEGER(KIND=JPIM),                   INTENT(IN)    :: KBLOCK   ! specify a block for threaded access
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(OUT)   :: KBLOCK_SIZE ! Number of horizontal points (our work allocation)

INTEGER(KIND=JPIM) :: IFIELDS_USED, ID, IUSED_FIDS(SELF%NFIELDS),J
LOGICAL :: LLON,LLHDERS
!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------
! Do not add DR HOOK to this routine: it causes poor performance
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_COMMON',0,ZHOOK_HANDLE)

ACCESS_CONTROL%IBLOCK = KBLOCK

! Set up a list of the fields in the access type, depending on what
! was requested (fids) and whether those fields are available in the
! fields container.
IFIELDS_USED = 0
IF (PRESENT(KFIDS)) THEN
  DO J=1,SIZE(KFIDS)
    ID = KFIDS(J)
    LLON = .FALSE.
    LLHDERS = .FALSE.
    IF (SELF%FIELD_INDEX(ID)%L_ON) THEN
      LLON = .TRUE.
      IF(SELF%INDEX_GP(ID)%L_HDERS) LLHDERS = .TRUE.
    ENDIF
    IF(LLON) THEN
      IFIELDS_USED = IFIELDS_USED + 1
      IUSED_FIDS(IFIELDS_USED) = ID
    ENDIF

  ENDDO
ELSE

  DO ID=1,SELF%NFIELDS

    LLON = .FALSE.
    LLHDERS = .FALSE.
    IF (SELF%FIELD_INDEX(ID)%L_ON) THEN
      LLON = .TRUE.
      IF(SELF%INDEX_GP(ID)%L_HDERS) LLHDERS = .TRUE.
    ENDIF

    IF(LLON) THEN
      IFIELDS_USED = IFIELDS_USED + 1
      IUSED_FIDS(IFIELDS_USED) = ID
    ENDIF
  ENDDO
ENDIF

IF(PRESENT(KFIDS)) THEN
  IF(IFIELDS_USED < SIZE(KFIDS)) THEN
    WRITE(0,*) 'FIELD_OPEN - KFIDS MISSING '
    DO J=1,SIZE(KFIDS)
      IF(.NOT. ANY(IUSED_FIDS(1:IFIELDS_USED) == KFIDS(J))) WRITE(0,*) 'FID ',J,KFIDS(J),' MISSING'
    ENDDO
  ENDIF
ENDIF


ALLOCATE(KFIELDS(IFIELDS_USED))
KFIELDS(:) = IUSED_FIDS(1:IFIELDS_USED)

IF(PRESENT(KBLOCK_SIZE)) KBLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_COMMON',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_COMMON

!=======================================================================================
SUBROUTINE FIELD_CLOSE_FAC(SELF, FAC)
! Finish accessing model fields (access type):
!  - Release the "lock" on these fields
!  - Nullify all the access pointers, so they cannot
! be used any more (at least, not without a seg fault!)
! -------------------------------------------------
CLASS(FIELD_CONTAINER_GP), INTENT(INOUT) :: SELF
CLASS(FIELD_ACCESS_BASE) , INTENT(INOUT) :: FAC
!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Do not add DR HOOK to this routine: it causes poor performance
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_FAC',0,ZHOOK_HANDLE)

CALL MAP_NAMES_TO_STORAGE(SELF, FAC%FIELDS, FAC%ACCESS_CONTROL%IBLOCK, MFAC=FAC, LD_RELEASE=.TRUE.)

SELECT TYPE (FAC)
TYPE IS (FIELD_ACCESS)
  IF(SELF%LHDERS) THEN
    IF(ASSOCIATED(FAC%DL)) DEALLOCATE(FAC%DL)
    IF(ASSOCIATED(FAC%DM)) DEALLOCATE(FAC%DM)
  ENDIF
  DEALLOCATE(FAC%FIELDS)
END SELECT
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_FAC',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_FAC
!=======================================================================================
SUBROUTINE FIELD_OPEN_ONE_1D(SELF, KFID, KBLOCK, PFIELD, KBLOCK_SIZE)
! -------------------------------------------------
! Get access to one model field via a pointer
! There are versions for 1D and 2D fields
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),  INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),         INTENT(IN)    :: KFID      ! field ID
INTEGER(KIND=JPIM),         INTENT(IN)    :: KBLOCK   ! specify a block for threaded access
REAL(KIND=JPRB), POINTER,   INTENT(OUT)   :: PFIELD(:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KBLOCK_SIZE ! Number of horizontal points (our work allocation)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_ONE_1D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

CALL FMAP_1D(SELF, KBLOCK, KFID, PFIELD)

IF(PRESENT(KBLOCK_SIZE)) KBLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_ONE_1D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_ONE_1D
!=======================================================================================
SUBROUTINE FIELD_OPEN_ONE_2D(SELF, KFID, KBLOCK, PFIELD, KBLOCK_SIZE)

CLASS(FIELD_CONTAINER_GP),    INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),           INTENT(IN)    :: KFID ! field ID
INTEGER(KIND=JPIM),           INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB), POINTER,     INTENT(OUT)   :: PFIELD(:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT)   :: KBLOCK_SIZE ! Number of horizontal points (our work allocation)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_ONE_2D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

CALL FMAP_2D(SELF, KBLOCK, KFID, PFIELD)

IF(PRESENT(KBLOCK_SIZE)) KBLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_ONE_2D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_ONE_2D
!=======================================================================================
SUBROUTINE FIELD_OPEN_ONE_3D(SELF, KFID, KBLOCK, PFIELD, KBLOCK_SIZE)

CLASS(FIELD_CONTAINER_GP),    INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),           INTENT(IN)    :: KFID ! field ID
INTEGER(KIND=JPIM),           INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB), POINTER,     INTENT(OUT)   :: PFIELD(:,:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT)   :: KBLOCK_SIZE ! Number of horizontal points (our work allocation)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_ONE_3D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

CALL FMAP_3D(SELF, KBLOCK, KFID, PFIELD)

IF(PRESENT(KBLOCK_SIZE)) KBLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_ONE_3D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_ONE_3D
!=======================================================================================
SUBROUTINE FIELD_CLOSE_ONE_1D(SELF, KFID, KBLOCK, PFIELD)
! -------------------------------------------------
! Finish accessing model fields (one field)
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)    :: KFID ! field ID
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK   ! specify a block for threaded access
REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: PFIELD(:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_ONE_1D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

NULLIFY(PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_ONE_1D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_ONE_1D
!=======================================================================================
SUBROUTINE FIELD_CLOSE_ONE_2D(SELF, KFID, KBLOCK, PFIELD)

CLASS(FIELD_CONTAINER_GP),INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)    :: KFID ! field ID
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: PFIELD(:,:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_ONE_2D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

IF( SELF%MSTORTYPE == JP_STORAGEOLD .AND. KBLOCK == SELF%NBLOCKS .AND. &
              &  .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN
  ! Special last block hack for GFL interoperability only.
  ! (we need a copy-out and copy-in approach)
  CALL UNMAP_LAST_BLOCK(SELF, KBLOCK, KFID)
ENDIF

NULLIFY(PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_ONE_2D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_ONE_2D
!=======================================================================================
SUBROUTINE FIELD_CLOSE_ONE_3D(SELF, KFID, KBLOCK, PFIELD)

CLASS(FIELD_CONTAINER_GP),INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)    :: KFID ! field ID
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: PFIELD(:,:,:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_ONE_3D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

IF( SELF%MSTORTYPE == JP_STORAGEOLD .AND. KBLOCK == SELF%NBLOCKS .AND. &
              &  .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN
  ! Special last block hack for GFL interoperability only.
  ! (we need a copy-out and copy-in approach)
  CALL UNMAP_LAST_BLOCK(SELF, KBLOCK, KFID)
ENDIF

NULLIFY(PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_ONE_3D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_ONE_3D
!=======================================================================================
SUBROUTINE FIELD_OPEN_MULTI(SELF, KFIDS, KBLOCK, PFIELD, KBLOCK_SIZE)
! --------------------------------------------------
! Get access to a group of 2D model fields as a 3D block
! --------------------------------------------------

CLASS(FIELD_CONTAINER_GP),    INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),           INTENT(IN)    :: KFIDS(:) ! FIELD ID
INTEGER(KIND=JPIM),           INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB), POINTER,     INTENT(OUT)   :: PFIELD(:,:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT)   :: KBLOCK_SIZE ! Number of horizontal points (our work allocation)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_MULTI',0,ZHOOK_HANDLE)
CALL FMAP_MULTI(SELF, KBLOCK, KFIDS, PFIELD)

IF(PRESENT(KBLOCK_SIZE)) THEN
  KBLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)
ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_OPEN_MULTI',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_MULTI
!=======================================================================================

SUBROUTINE FIELD_CLOSE_MULTI(SELF, KFIDS, KBLOCK, PFIELD)
! -------------------------------------------------
! Finish accessing model fields
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)    :: KFIDS(:) ! field ID
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: PFIELD(:,:,:)

INTEGER(KIND=JPIM) :: ISIZE_GFL, IPOINTS, ILEVELS, IFIELDS
REAL(KIND=JPRB), POINTER :: ZTMP(:,:,:) => NULL()
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_MULTI',0,ZHOOK_HANDLE)
IF( SELF%MSTORTYPE == JP_STORAGEOLD .AND. KBLOCK == SELF%NBLOCKS .AND. .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

  ! Special last block hack for GFL interoperability only.
  ! (we need a copy-out and copy-in approach)
  CALL GFL_AS_A_WHOLE(KBLOCK, SELF%M_ATTACH,SELF%YRGFL, ZTMP)
  IPOINTS = SELF%NPOINTS_BLOCK(KBLOCK)
  ILEVELS = SELF%NLEVELS(SELF%METADATA(KFIDS(1))%D2TYPE) ! NFTBC assume all fids same size
  IFIELDS = MIN(SIZE(SELF%STORAGE_LAST_BLOCK,3),SIZE(ZTMP,3)) ! AJGDB really really hacky now
  ZTMP(1:IPOINTS,1:ILEVELS,1:IFIELDS) = SELF%STORAGE_LAST_BLOCK(:,1:ILEVELS,1:IFIELDS)

#endif

ENDIF

NULLIFY(PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_CLOSE_MULTI',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_MULTI


FUNCTION FIELD_ITERATED(SELF, KFIDS, KBLOCK, FAC,LDHDERS,KOTID,KBLOCK_SIZE,FVAR)
! -------------------------------------------------
! Iterator (do while) access to model fields.
!
! - Given a list of field IDs, each successive call
!   to this function will return the next field.
!
! - All fields in the list should have the same shape
!   (e.g. all 1D, all 2D)
!
! Function return value is 0 if no more fields are
! available and 1 otherwise.
!
! The reason for choosing a do while loop is it allows
! us to control the pointer access (i.e. to nullify
! pointers into the field storage so they are not
! used outside the scope of the do while loop)

INTEGER(KIND=JPIM) :: FIELD_ITERATED

CLASS(FIELD_CONTAINER_GP),TARGET,     INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, TARGET, INTENT(IN)    :: KFIDS(:)      ! list of field IDs
INTEGER(KIND=JPIM),                   INTENT(IN)    :: KBLOCK        ! block specifier for threaded access
CLASS(FIELD_ACCESS_BASE),             INTENT(INOUT) :: FAC
LOGICAL,OPTIONAL,                     INTENT(IN)    :: LDHDERS
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(OUT)   :: KOTID         ! ID of the field being returned
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(OUT)   :: KBLOCK_SIZE   ! Number of horizontal points (our work allocation)
TYPE(TYPE_FVAR),POINTER,OPTIONAL,     INTENT(OUT)   :: FVAR
INTEGER(KIND=JPIM)          :: ID, INO_FLDS_ITER
INTEGER(KIND=JPIM), POINTER :: IUSE_FIDS(:)
LOGICAL :: LLHDERS
!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_ITERATED',0,ZHOOK_HANDLE)

IF(PRESENT(LDHDERS)) THEN
  LLHDERS = LDHDERS
ELSE
  LLHDERS = .FALSE.
ENDIF
SELECT TYPE (FAC)
TYPE IS (FIELD_ACCESS)
  IF(LLHDERS) THEN
    IF(.NOT. ASSOCIATED(FAC%DL)) ALLOCATE(FIELD_ACCESS::FAC%DL)
    IF(.NOT. ASSOCIATED(FAC%DM)) ALLOCATE(FIELD_ACCESS::FAC%DM)
  ENDIF
END SELECT

IF(PRESENT(KFIDS)) THEN
  IUSE_FIDS=>KFIDS ! AJGDB - is this guaranteed not to disapper? No!
ELSE
  IUSE_FIDS=>SELF%MALL_FIELDS
ENDIF
INO_FLDS_ITER = SIZE(IUSE_FIDS)

IF(FAC%ITERATOR_COUNT < INO_FLDS_ITER) THEN

  FAC%ITERATOR_COUNT = FAC%ITERATOR_COUNT+1
  ID = IUSE_FIDS(FAC%ITERATOR_COUNT)
  IF(PRESENT(KOTID)) KOTID = ID
  IF(PRESENT(FVAR)) FVAR=>SELF%METADATA(ID)

  IF(SELF%FIELD_INDEX(ID)%L_ON) THEN
    ! Give access to one field via a pointer
    IF(SELF%METADATA(ID)%NDIMS==1) THEN
      IF(LLHDERS) THEN
        SELECT TYPE (FAC)
        TYPE IS (FIELD_ACCESS)
          CALL FMAP_1D(SELF, KBLOCK, ID, FAC%RR1,.TRUE., FAC%DL%RR1,FAC%DM%RR1)
        END SELECT
      ELSE
        CALL FMAP_1D(SELF, KBLOCK, ID, FAC%RR1)
      ENDIF
      NULLIFY(FAC%RR2)
    ELSEIF(SELF%METADATA(ID)%NDIMS==2) THEN
      IF(LLHDERS) THEN
        SELECT TYPE (FAC)
        TYPE IS (FIELD_ACCESS)
          CALL FMAP_2D(SELF, KBLOCK, ID, FAC%RR2,.TRUE., FAC%DL%RR2,FAC%DM%RR2)
        END SELECT
      ELSE
         CALL FMAP_2D(SELF, KBLOCK, ID, FAC%RR2)
       ENDIF
      NULLIFY(FAC%RR1)
    ENDIF
    FIELD_ITERATED = 1
  ELSE
    NULLIFY(FAC%RR2)
    NULLIFY(FAC%RR1)
    FIELD_ITERATED = 0
  ENDIF
ELSE
  NULLIFY(FAC%RR2)
  NULLIFY(FAC%RR1)
  FIELD_ITERATED = 0
  FAC%ITERATOR_COUNT = 0
ENDIF

IF(PRESENT(KBLOCK_SIZE)) KBLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_ITERATED',1,ZHOOK_HANDLE)

END FUNCTION FIELD_ITERATED
!=======================================================================================
SUBROUTINE FIELD_EQUALS(SELF, PVALUE, KFIDS, KBLOCK)
! -------------------------------------------------
! Set fields to a value
! -------------------------------------------------
CLASS(FIELD_CONTAINER_GP),TARGET,     INTENT(INOUT) :: SELF
REAL(KIND=JPRB),                      INTENT(IN)    :: PVALUE
INTEGER(KIND=JPIM), OPTIONAL, TARGET, INTENT(IN)    :: KFIDS(:)
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(IN)    :: KBLOCK

INTEGER(KIND=JPIM) :: I, ID, ISTART, IEND
INTEGER(KIND=JPIM) :: IUSE_BLOCK, IBLOCK, IBLOCK_START, IBLOCK_END
INTEGER(KIND=JPIM), POINTER :: IUSE_FIDS(:)
REAL(KIND=JPRB), POINTER :: F1D(:), F2D(:,:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_EQUALS',0,ZHOOK_HANDLE)
IF(PRESENT(KBLOCK)) THEN
  ! Access from within a threaded region
  IUSE_BLOCK = KBLOCK
  IBLOCK_START = KBLOCK
  IBLOCK_END   = KBLOCK
ELSE
  ! Non-threaded access
  IUSE_BLOCK = JP_ALL_BLOCKS
  IBLOCK_START = 1
  IBLOCK_END = SELF%NBLOCKS
ENDIF

IF(PRESENT(KFIDS)) THEN
  IUSE_FIDS => KFIDS
ELSE
  IUSE_FIDS => SELF%MALL_FIELDS
ENDIF


  ! NFTBC CAN BE OMP PARALLEL EVENTUALLY
DO IBLOCK = IBLOCK_START, IBLOCK_END
  DO I=1,SIZE(IUSE_FIDS)
    ID = IUSE_FIDS(I)
    IF(SELF%FIELD_INDEX(ID)%L_ON) THEN
      IF (SELF%METADATA(ID)%NDIMS == 1) THEN

        CALL FMAP_1D(SELF, IBLOCK, ID, F1D)
        F1D = PVALUE

      ELSEIF (SELF%METADATA(ID)%NDIMS == 2) THEN

        CALL FMAP_2D(SELF, IBLOCK, ID, F2D)
        F2D = PVALUE

        IF( SELF%MSTORTYPE == JP_STORAGEOLD .AND. IBLOCK == SELF%NBLOCKS .AND. &
         &  .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN
            ! SPECIAL LAST BLOCK HACK FOR GFL INTEROPERABILITY ONLY.
          CALL UNMAP_LAST_BLOCK(SELF, IBLOCK, ID)
        ENDIF

      ENDIF
    ENDIF
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_EQUALS',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_EQUALS
!=======================================================================================
SUBROUTINE FIELD_ADD(F1, F2, KFIDS, KBLOCK)
! -------------------------------------------------
! Add fields: f1=f1+f2
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),                INTENT(INOUT) :: F1
CLASS(FIELD_CONTAINER_GP),                INTENT(INOUT) :: F2
INTEGER(KIND=JPIM), OPTIONAL, TARGET, INTENT(IN)    :: KFIDS(:)
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(IN)    :: KBLOCK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_ADD',0,ZHOOK_HANDLE)
CALL TWO_FIELD_ACTION(F1, F2, JP_FADD, KFIDS, KBLOCK)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_ADD',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_ADD
!=======================================================================================

SUBROUTINE FIELD_COPY(F1, F2, KFIDS, KBLOCK)
! -------------------------------------------------
! Copy fields: f1=f2
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),                INTENT(INOUT) :: F1
CLASS(FIELD_CONTAINER_GP),                INTENT(INOUT) :: F2
INTEGER(KIND=JPIM), OPTIONAL, TARGET, INTENT(IN)    :: KFIDS(:)
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(IN)    :: KBLOCK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_COPY',0,ZHOOK_HANDLE)
CALL TWO_FIELD_ACTION(F1, F2, JP_FCOPY, KFIDS, KBLOCK)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_COPY',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_COPY
!=======================================================================================

SUBROUTINE TWO_FIELD_ACTION(F1, F2, KACTION, KFIDS, KBLOCK)
! -------------------------------------------------
! Two-field actions:
!
!    add:   f1 = f1+f2
!    copy:  f1 = f2
!
!  - f1 and f2 need not be identical: only the
!    intersection of fields in f1 and f2 is considered
!  - If a list of "fids" is supplied, only the fields
!    in this list are considered
!  - If "block" is supplied, only one block is considered
!
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),TARGET,         INTENT(INOUT) :: F1
CLASS(FIELD_CONTAINER_GP),                INTENT(INOUT) :: F2
INTEGER(KIND=JPIM),                   INTENT(IN)    :: KACTION
INTEGER(KIND=JPIM), OPTIONAL, TARGET, INTENT(IN)    :: KFIDS(:)
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(IN)    :: KBLOCK

INTEGER(KIND=JPIM) :: ID
INTEGER(KIND=JPIM) :: IUSE_BLOCK, IBLOCK, IBLOCK_START, IBLOCK_END
INTEGER(KIND=JPIM), POINTER :: IUSE_FIDS(:)
REAL(KIND=JPRB), POINTER :: F1_1D(:), F1_2D(:,:), F2_1D(:), F2_2D(:,:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:TWO_FIELD_ACTION',0,ZHOOK_HANDLE)
IF(PRESENT(KBLOCK)) THEN
  ! ACCESS FROM WITHIN A THREADED REGION
  IUSE_BLOCK = KBLOCK
  IBLOCK_START = KBLOCK
  IBLOCK_END   = KBLOCK
ELSE
  ! NON-THREADED ACCESS
  IUSE_BLOCK = JP_ALL_BLOCKS
  IBLOCK_START = 1
  IBLOCK_END = F1%NBLOCKS
ENDIF

IF(PRESENT(KFIDS)) THEN
  IUSE_FIDS => KFIDS
ELSE
  IUSE_FIDS => F1%MALL_FIELDS
ENDIF


! NFTBC CAN BE OMP PARALLEL EVENTUALLY
DO IBLOCK = IBLOCK_START, IBLOCK_END
  DO ID = 1, F1%NFIELDS
    IF (F1%FIELD_INDEX(ID)%L_ACTIVE .AND. F2%FIELD_INDEX(ID)%L_ACTIVE .AND. ANY(ID == IUSE_FIDS)) THEN
      IF (F1%METADATA(ID)%NDIMS == 1) THEN

        CALL FMAP_1D(F1, IBLOCK, ID, F1_1D)
        CALL FMAP_1D(F2, IBLOCK, ID, F2_1D)
        IF     (KACTION==JP_FADD)  THEN
          F1_1D = F1_1D + F2_1D
        ELSEIF (KACTION==JP_FCOPY) THEN
          F1_1D = F2_1D
        ENDIF

      ELSEIF (F1%METADATA(ID)%NDIMS == 2) THEN

        CALL FMAP_2D(F1, IBLOCK, ID, F1_2D)
        CALL FMAP_2D(F2, IBLOCK, ID, F2_2D)
        IF     (KACTION==JP_FADD)  THEN
          F1_2D = F1_2D + F2_2D
        ELSEIF (KACTION==JP_FCOPY) THEN
          F1_2D = F2_2D
        ENDIF

        IF( F1%MSTORTYPE == JP_STORAGEOLD .AND. IBLOCK == F1%NBLOCKS .AND. &
         &  .NOT. F1%L_ATTACH_GFL_NPROMA_FIXED) THEN
          ! SPECIAL LAST BLOCK HACK FOR GFL INTEROPERABILITY ONLY.
          CALL UNMAP_LAST_BLOCK(F1, IBLOCK, ID)
        ENDIF

      ENDIF
    ENDIF
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:TWO_FIELD_ACTION',1,ZHOOK_HANDLE)


END SUBROUTINE TWO_FIELD_ACTION
!=======================================================================================
SUBROUTINE FIELD_NORM(SELF)
! -------------------------------------------------
! Compute and print norms
! -------------------------------------------------
CLASS(FIELD_CONTAINER_GP), INTENT(INOUT) :: SELF

INTEGER(KIND=JPIM) :: ID, ISIZE, IBLOCK
REAL(KIND=JPRB), POINTER :: ZVAR1D(:), ZVAR2D(:,:), ZVAR3D(:,:,:)
REAL(KIND=JPRB) :: ZRMEAN, ZRRMS, ZRMIN, ZRMAX, ZMEAN, ZRMS, ZMIN, ZMAX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_NORM',0,ZHOOK_HANDLE)
WRITE(NULOUT,*) '  Name    Ndims     Mean         RMS         Max        Min'
DO ID=1,SELF%NFIELDS
  IF (SELF%FIELD_INDEX(ID)%L_ON) THEN

    ! Norms are per-PE for the moment. Also, the algorithm is simplistic
    ! (e.g. prone to overflow?) right now. NFTBC
    ZRMEAN=0.0_JPRB
    ZRRMS =0.0_JPRB
    ZRMIN =HUGE(ZRMIN)
    ZRMAX =-HUGE(ZRMAX)
    DO IBLOCK = 1, SELF%NBLOCKS
      IF (SELF%METADATA(ID)%NDIMS == 1) THEN

        CALL FMAP_1D(SELF, IBLOCK, ID, ZVAR1D)
        ZMEAN= SUM(ZVAR1D)
        ZRMS = SUM(ZVAR1D**2)
        ZMIN = MINVAL(ZVAR1D)
        ZMAX = MAXVAL(ZVAR1D)

      ELSEIF  (SELF%METADATA(ID)%NDIMS == 2) THEN

        CALL FMAP_2D(SELF, IBLOCK, ID, ZVAR2D)
        ZMEAN= SUM(ZVAR2D)
        ZRMS = SUM(ZVAR2D**2)
        ZMIN = MINVAL(ZVAR2D)
        ZMAX = MAXVAL(ZVAR2D)

      ELSEIF  (SELF%METADATA(ID)%NDIMS == 3) THEN

        CALL FMAP_3D(SELF, IBLOCK, ID, ZVAR3D)
        ZMEAN= SUM(ZVAR3D)
        ZRMS = SUM(ZVAR3D**2)
        ZMIN = MINVAL(ZVAR3D)
        ZMAX = MAXVAL(ZVAR3D)

      ENDIF
      ZRMEAN = ZRMEAN + ZMEAN
      ZRRMS  = ZRRMS  + ZRMS
      ZRMIN  = MIN(ZRMIN,ZMIN)
      ZRMAX  = MAX(ZRMAX,ZMAX)
    ENDDO

    IF (SELF%METADATA(ID)%NDIMS == 1) THEN
      ISIZE = SELF%NPOINTS
    ELSE
      ISIZE = SELF%NPOINTS*SELF%NLEVELS(SELF%METADATA(ID)%D2TYPE)
    ENDIF

    ZRMEAN = ZRMEAN / ISIZE
    ZRRMS  = SQRT(ZRRMS / ISIZE)

    WRITE(NULOUT,"(A10,X,I3,X,4F12.5)") &
      & SELF%METADATA(ID)%NAME, SELF%METADATA(ID)%NDIMS, ZRMEAN, ZRRMS, ZRMAX, ZRMIN

  ENDIF
ENDDO
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FIELD_NORM',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_NORM
!=======================================================================================
FUNCTION FIELD_NBLOCKS(SELF)
! ------------------------------------------------------
! Find the total number of blocks
! ------------------------------------------------------
INTEGER(KIND=JPIM)                :: FIELD_NBLOCKS

CLASS(FIELD_CONTAINER_GP), INTENT(IN) :: SELF

FIELD_NBLOCKS = SELF%NBLOCKS

END FUNCTION FIELD_NBLOCKS
!=======================================================================================
FUNCTION FIELD_BLOCK_SIZE(SELF,KBLOCK)
! ------------------------------------------------------
! Find the horizontal dimension of a particular block
! ------------------------------------------------------

INTEGER(KIND=JPIM)                :: FIELD_BLOCK_SIZE
CLASS(FIELD_CONTAINER_GP), INTENT(IN) :: SELF
INTEGER(KIND=JPIM),        INTENT(IN) :: KBLOCK

FIELD_BLOCK_SIZE = SELF%NPOINTS_BLOCK(KBLOCK)

END FUNCTION FIELD_BLOCK_SIZE
!=======================================================================================

SUBROUTINE MAP_NAMES_TO_STORAGE(SELF, KFIDS, KBLOCK, MFAC, LD_RELEASE,KFIDS_HD)
! -------------------------------------------------
! Fill up the field_access type
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),     INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),            INTENT(IN)    :: KFIDS(:), KBLOCK
CLASS(FIELD_ACCESS_BASE), OPTIONAL, INTENT(INOUT) :: MFAC
LOGICAL, OPTIONAL,             INTENT(IN)    :: LD_RELEASE ! Indicates when access is being closed
INTEGER(KIND=JPIM), TARGET, OPTIONAL,INTENT(IN) :: KFIDS_HD(:)

INTEGER(KIND=JPIM)       :: ID,INFIDS,INHFIDS,IFIDS(SELF%NFIELDS),J
LOGICAL                  :: LL_NULLIFY, LL_RELEASE,LLHDERS,LLHFIDS
REAL(KIND=JPRB), POINTER :: FIELD_1D(:), FIELD_1D_DL(:),FIELD_1D_DM(:)
REAL(KIND=JPRB), POINTER :: FIELD_2D(:,:),FIELD_2D_DL(:,:),FIELD_2D_DM(:,:)
REAL(KIND=JPRB), POINTER :: FIELD_3D(:,:,:),FIELD_3D_DL(:,:,:),FIELD_3D_DM(:,:,:)
!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Do not add DR HOOK to this routine: it causes poor performance
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:MAP_NAMES_TO_STORAGE',0,ZHOOK_HANDLE)
LL_RELEASE = .FALSE.
IF(PRESENT(LD_RELEASE)) LL_RELEASE = LD_RELEASE
LLHFIDS = PRESENT(KFIDS_HD)
INFIDS = SIZE(KFIDS)
IFIDS(1:INFIDS) = KFIDS(:)
IF(LLHFIDS) THEN
  INHFIDS=SIZE(KFIDS_HD)
  IFIDS(INFIDS+1:INFIDS+INHFIDS) = KFIDS_HD(:)
ELSE
  INHFIDS=0
ENDIF

DO J=1,INFIDS+INHFIDS
  ID = IFIDS(J)
  IF(LLHFIDS) THEN
    LLHDERS = SELF%INDEX_GP(ID)%L_HDERS .AND. J>INFIDS
  ELSE
    LLHDERS = .FALSE.
  ENDIF
  LL_NULLIFY = .FALSE.

  IF( .NOT. LL_RELEASE) THEN

    IF (SELF%METADATA(ID)%NDIMS == 1) THEN
      CALL FMAP_1D(SELF, KBLOCK, ID, FIELD_1D,LLHDERS, FIELD_1D_DL,FIELD_1D_DM)
    ELSEIF (SELF%METADATA(ID)%NDIMS == 2) THEN
      CALL FMAP_2D(SELF, KBLOCK, ID, FIELD_2D,LLHDERS, FIELD_2D_DL,FIELD_2D_DM)
    ELSEIF (SELF%METADATA(ID)%NDIMS == 3) THEN
      CALL FMAP_3D(SELF, KBLOCK, ID, FIELD_3D,LLHDERS, FIELD_3D_DL,FIELD_3D_DM)
    ELSE
      CALL ABOR1('FIELD_CONTAINER_GP_MOD: Unsupported number of dimensions')
    ENDIF

  ELSEIF( LL_RELEASE) THEN

    ! "Unmapping" is only a temporary thing while we have to maintain
    ! inter-operability with the GFLs etc. And only for the last block of
    ! 2d fields. NFTBC
    IF( SELF%MSTORTYPE == JP_STORAGEOLD .AND. SELF%METADATA(ID)%NDIMS == 2 .AND. KBLOCK == SELF%NBLOCKS .AND. &
              &  .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN
      CALL UNMAP_LAST_BLOCK(SELF, KBLOCK, ID)
    ENDIF

    LL_NULLIFY = .TRUE.

  ELSE
    ! NFTBC /AJGDB we can probably avoid doing this for unused fields: they should have
    ! started nullified (from the type definition => null()) and stayed nullified so
    ! this is completely unnecessary.
    LL_NULLIFY = .TRUE.
  ENDIF
  IF(.NOT.LL_NULLIFY) THEN

    CALL MFAC%FIELD_MAP_STORAGE( ID, STORAGE_1D=FIELD_1D, STORAGE_2D=FIELD_2D, STORAGE_3D=FIELD_3D, &
     & LD_NULLIFY=LL_NULLIFY)
    SELECT TYPE (MFAC)
      TYPE IS (FIELD_ACCESS)
      IF(LLHDERS) THEN
        CALL MFAC%DL%FIELD_MAP_STORAGE( ID, STORAGE_1D=FIELD_1D_DL, STORAGE_2D=FIELD_2D_DL, STORAGE_3D=FIELD_3D_DL, &
         & LD_NULLIFY=LL_NULLIFY)
        CALL MFAC%DM%FIELD_MAP_STORAGE( ID, STORAGE_1D=FIELD_1D_DM, STORAGE_2D=FIELD_2D_DM, STORAGE_3D=FIELD_3D_DM, &
         & LD_NULLIFY=LL_NULLIFY)
      ENDIF
    END SELECT

  ENDIF
ENDDO

!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:MAP_NAMES_TO_STORAGE',1,ZHOOK_HANDLE)

END SUBROUTINE MAP_NAMES_TO_STORAGE
!=======================================================================================
SUBROUTINE FMAP_1D(SELF, KBLOCK, KID, FPTR,LDHDERS,FPTR_L,FPTR_M)
! -------------------------------------------------
! Maps a pointer into the storage area (1D)
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),    INTENT(IN),TARGET  :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)  :: KBLOCK, KID
REAL(KIND=JPRB), POINTER, INTENT(OUT) :: FPTR(:)
LOGICAL,  OPTIONAL,     INTENT(IN)    :: LDHDERS
REAL(KIND=JPRB), OPTIONAL,POINTER,INTENT(OUT) :: FPTR_L(:)
REAL(KIND=JPRB), OPTIONAL,POINTER,INTENT(OUT) :: FPTR_M(:)

INTEGER(KIND=JPIM) :: ISTART, IEND, IPOINTS
LOGICAL :: LL_WANTED,LLHDERS
!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_1D',0,ZHOOK_HANDLE)

LLHDERS = .FALSE.
IF(PRESENT(LDHDERS)) LLHDERS = LDHDERS
IF(SELF%METADATA(KID)%NDIMS == 1) THEN
  IF(SELF%FIELD_INDEX(KID)%L_ON) THEN

    IPOINTS = SELF%NPOINTS_BLOCK(KBLOCK)

    IF (SELF%MSTORTYPE == JP_STORAGE1D) THEN

      ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%ISTART
      IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%IEND

      FPTR => SELF%STORAGE1D(ISTART:IEND)
      IF(LLHDERS) THEN
        IF(PRESENT(FPTR_L)) THEN
          ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%DL%ISTART
          IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%DL%IEND
          FPTR_L => SELF%STORAGE1D(ISTART:IEND)
        ENDIF
        IF(PRESENT(FPTR_M)) THEN
          ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%DM%ISTART
          IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%DM%IEND
          FPTR_M => SELF%STORAGE1D(ISTART:IEND)
        ENDIF
      ENDIF
    ELSEIF(SELF%MSTORTYPE == JP_STORAGE_ARR) THEN
      FPTR => SELF%STORAGE_ARRAYS(KID)%STORE2D(:,KBLOCK)
    ELSEIF (SELF%MSTORTYPE == JP_STORAGE3D) THEN

      !AJGDB

    ELSEIF (SELF%MSTORTYPE == JP_STORAGEOLD) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

      IF(SELF%L_ATTACH_GFL_NPROMA_FIXED)IPOINTS=MAXVAL(SELF%NPOINTS_BLOCK)

      CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, 1, KID, KBLOCK, SELF%M_ATTACH, &
      & SELF%YRGMV,SELF%YRGFL,SELF%YRSURF,FPTR_1D=FPTR, &
      & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
      IF(LLHDERS) THEN
        IF(SELF%M_ATTACH == JP_T0) THEN
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, 1, KID, KBLOCK, JP_T0_DL, &
           & SELF%YRGMV,SELF%YRGFL,SELF%YRSURF,FPTR_1D=FPTR_L, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, 1, KID, KBLOCK, JP_T0_DM, &
           & SELF%YRGMV,SELF%YRGFL,SELF%YRSURF,FPTR_1D=FPTR_M, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
        ELSEIF(SELF%M_ATTACH == JP_T9) THEN
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, 1, KID, KBLOCK, JP_T9_DL, &
           & SELF%YRGMV,SELF%YRGFL,SELF%YRSURF,FPTR_1D=FPTR_L, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, 1, KID, KBLOCK, JP_T9_DM, &
           & SELF%YRGMV,SELF%YRGFL,SELF%YRSURF,FPTR_1D=FPTR_M, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
        ELSE
          CALL ABOR1('FIELD_CONTAINER_GP_MOD:FMAP_1D - NO DERIVATIVES')
        ENDIF
      ENDIF

#endif

    ENDIF

  ELSE
    NULLIFY(FPTR)
  ENDIF
ELSE
  CALL ABOR1('FIELD_CONTAINER_GP_MOD: '//SELF%METADATA(KID)%NAME//' IS NOT 1D')
  NULLIFY(FPTR)
ENDIF
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_1D',1,ZHOOK_HANDLE)

END SUBROUTINE FMAP_1D
!=======================================================================================

SUBROUTINE FMAP_2D(SELF, KBLOCK, KID, FPTR,LDHDERS,FPTR_L,FPTR_M)
! -------------------------------------------------
! Maps a pointer into the storage area (2D)
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),TARGET, INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),               INTENT(IN)    :: KBLOCK, KID
REAL(KIND=JPRB),          POINTER,INTENT(OUT)   :: FPTR(:,:)
LOGICAL,         OPTIONAL,        INTENT(IN)    :: LDHDERS
REAL(KIND=JPRB), OPTIONAL,POINTER,INTENT(OUT)   :: FPTR_L(:,:)
REAL(KIND=JPRB), OPTIONAL,POINTER,INTENT(OUT)   :: FPTR_M(:,:)

INTEGER(KIND=JPIM) :: ISTART, IEND, IPOINTS, I2D, ILEVELS
REAL(KIND=JPRB), POINTER :: ZNEC_WORKAROUND(:) ! HELPS AVOID A COMPILER BUG ON THE NEC
LOGICAL :: LL_WANTED,LLHDERS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_2D',0,ZHOOK_HANDLE)

LLHDERS = .FALSE.
IF(PRESENT(LDHDERS)) LLHDERS = LDHDERS
IF(SELF%METADATA(KID)%NDIMS == 2) THEN
  IF(SELF%FIELD_INDEX(KID)%L_ON) THEN

    IPOINTS = SELF%NPOINTS_BLOCK(KBLOCK)
    ILEVELS = SELF%NLEVELS(SELF%METADATA(KID)%D2TYPE)

    IF (SELF%MSTORTYPE == JP_STORAGE1D) THEN

      ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%ISTART
      IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%IEND

      FPTR(1:IPOINTS,1:ILEVELS) => SELF%STORAGE1D(ISTART:IEND)

      IF(LLHDERS) THEN
        IF(PRESENT(FPTR_L)) THEN
          ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%DL%ISTART
          IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%DL%IEND

          FPTR_L(1:IPOINTS,1:ILEVELS) => SELF%STORAGE1D(ISTART:IEND)

        ENDIF
        IF(PRESENT(FPTR_M)) THEN
          ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%DM%ISTART
          IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%DM%IEND

          FPTR_M(1:IPOINTS,1:ILEVELS) => SELF%STORAGE1D(ISTART:IEND)
        ENDIF
      ENDIF

    ELSEIF(SELF%MSTORTYPE == JP_STORAGE_ARR) THEN
      FPTR => SELF%STORAGE_ARRAYS(KID)%STORE3D(:,:,KBLOCK)
    ELSEIF (SELF%MSTORTYPE == JP_STORAGE3D) THEN

      !AJGDB

    ELSEIF (SELF%MSTORTYPE == JP_STORAGEOLD) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

      ! Direct pointer into the GFL etc...
      IF(SELF%L_ATTACH_GFL_NPROMA_FIXED) IPOINTS = MAXVAL(SELF%NPOINTS_BLOCK)
      CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, SELF%M_ATTACH, SELF%YRGMV, &
       & SELF%YRGFL,SELF%YRSURF,FPTR_2D=FPTR, FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
      IF(LLHDERS) THEN
        IF(SELF%M_ATTACH == JP_T0) THEN
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T0_DL, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_2D=FPTR_L, FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T0_DM, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_2D=FPTR_M, FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
        ELSEIF(SELF%M_ATTACH == JP_T9) THEN
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T9_DL, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_2D=FPTR_L, FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T9_DM, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_2D=FPTR_M, FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
        ELSE
          CALL ABOR1('field_container_gp_mod:fmap_2d - no derivatives')
        ENDIF
      ENDIF
      IF(KBLOCK == SELF%NBLOCKS .AND. ASSOCIATED(FPTR) .AND. .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN

        ! Special trick for the last, smaller, nproma block: take a copy of the
        ! correct size. Necessary because of the non rank-compliant F77 interfaces
        ! in the IFS (e.g. gprcp etc.). We can't give the user a strided 2D pointer.
        I2D = SELF%INDEX_LAST_BLOCK(KID)
        SELF%STORAGE_LAST_BLOCK(:,1:ILEVELS,I2D) = FPTR
        FPTR => SELF%STORAGE_LAST_BLOCK(:,1:ILEVELS,I2D)

      ENDIF

#endif

    ENDIF
  ELSE
    NULLIFY(FPTR)
  ENDIF
ELSE
  WRITE(0,*) 'FIELD_CONTAINER_GP_MOD: ',KID,' is not 2D'
  CALL ABOR1('FIELD_CONTAINER_GP_MOD: '//self%METADATA(kid)%name//' is not 2D')
  NULLIFY(FPTR)
ENDIF
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_2D',1,ZHOOK_HANDLE)

END SUBROUTINE FMAP_2D

SUBROUTINE FMAP_3D(SELF, KBLOCK, KID, FPTR,LDHDERS,FPTR_L,FPTR_M)
! -------------------------------------------------
! Maps a pointer into the storage area (3D)
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),TARGET, INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK, KID
REAL(KIND=JPRB), POINTER, INTENT(OUT)   :: FPTR(:,:,:)
LOGICAL,  OPTIONAL,     INTENT(IN)      :: LDHDERS
REAL(KIND=JPRB), OPTIONAL,POINTER,INTENT(OUT) :: FPTR_L(:,:,:)
REAL(KIND=JPRB), OPTIONAL,POINTER,INTENT(OUT) :: FPTR_M(:,:,:)


INTEGER(KIND=JPIM) :: ISTART, IEND, IPOINTS, I2D, ILEVELS,IDIM3
REAL(KIND=JPRB), POINTER :: ZNEC_WORKAROUND(:) ! HELPS AVOID A COMPILER BUG ON THE NEC
LOGICAL :: LL_WANTED,LLHDERS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_3D',0,ZHOOK_HANDLE)

LLHDERS = .FALSE.
IF(PRESENT(LDHDERS)) LLHDERS = LDHDERS
IF(SELF%METADATA(KID)%NDIMS == 3) THEN
  IF(SELF%FIELD_INDEX(KID)%L_ON) THEN

    IPOINTS = SELF%NPOINTS_BLOCK(KBLOCK)
    ILEVELS = SELF%NLEVELS(SELF%METADATA(KID)%D2TYPE)
    IDIM3   = SELF%NDIM3S(SELF%METADATA(KID)%D3TYPE)

    IF (SELF%MSTORTYPE == JP_STORAGE1D) THEN

      ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%ISTART
      IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%IEND

      FPTR(1:IPOINTS,1:ILEVELS,1:IDIM3) => SELF%STORAGE1D(ISTART:IEND)

      IF(LLHDERS) THEN
        IF(PRESENT(FPTR_L)) THEN
          ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%DL%ISTART
          IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%DL%IEND

          FPTR_L(1:IPOINTS,1:ILEVELS,1:IDIM3) => SELF%STORAGE1D(ISTART:IEND)

        ENDIF
        IF(PRESENT(FPTR_M)) THEN
          ISTART  = SELF%INDEX_BLOCK(KID,KBLOCK)%DM%ISTART
          IEND    = SELF%INDEX_BLOCK(KID,KBLOCK)%DM%IEND

          FPTR_M(1:IPOINTS,1:ILEVELS,1:IDIM3) => SELF%STORAGE1D(ISTART:IEND)
        ENDIF
      ENDIF

    ELSEIF(SELF%MSTORTYPE == JP_STORAGE_ARR) THEN
      FPTR =>  SELF%STORAGE_ARRAYS(KID)%STORE4D(:,:,:,KBLOCK)
    ELSEIF (SELF%MSTORTYPE == JP_STORAGE3D) THEN

      !AJGDB

    ELSEIF (SELF%MSTORTYPE == JP_STORAGEOLD) THEN

#ifndef FIELD_CONTAINER_GP_MOD_TEST

      ! Direct pointer into the GFL etc...
      IF(SELF%L_ATTACH_GFL_NPROMA_FIXED) IPOINTS = MAXVAL(SELF%NPOINTS_BLOCK)
      CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, SELF%M_ATTACH, SELF%YRGMV, &
                     & SELF%YRGFL,SELF%YRSURF,FPTR_3D=FPTR, &
                     & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL,KDIM3=IDIM3)
      IF(LLHDERS) THEN
        IF(SELF%M_ATTACH == JP_T0) THEN
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T0_DL, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_3D=FPTR_L, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL,KDIM3=IDIM3)
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T0_DM, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_3D=FPTR_M, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL,KDIM3=IDIM3)
        ELSEIF(SELF%M_ATTACH == JP_T9) THEN
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T9_DL, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_3D=FPTR_L, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL,KDIM3=IDIM3)
          CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, JP_T9_DM, SELF%YRGMV, &
           & SELF%YRGFL,SELF%YRSURF,FPTR_3D=FPTR_M, &
           & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL,KDIM3=IDIM3)
        ELSE
          CALL ABOR1('field_container_gp_mod:fmap_3d - no derivatives')
        ENDIF
      ENDIF

#endif

    ENDIF
  ELSE
    NULLIFY(FPTR)
  ENDIF
ELSE
  CALL ABOR1('FIELD_CONTAINER_GP_MOD: '//self%METADATA(kid)%name//' is not 3D')
  NULLIFY(FPTR)
ENDIF
!IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_3D',1,ZHOOK_HANDLE)

END SUBROUTINE FMAP_3D
!=======================================================================================
SUBROUTINE FMAP_MULTI(SELF, KBLOCK, KIDS, FPTR)
! -------------------------------------------------
! Maps a pointer into the storage area (3D -
! e.g. multiple 2D variables)
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),TARGET,INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK, KIDS(:)
REAL(KIND=JPRB), POINTER, INTENT(OUT)   :: FPTR(:,:,:)

INTEGER(KIND=JPIM) :: ISTART, IEND, IPOINTS, I2D, ILEVELS, ID, I, IVARS, ISIZE_GFL
REAL(KIND=JPRB), POINTER :: ZNEC_WORKAROUND(:) ! helps avoid a compiler bug on the NEC
LOGICAL :: LL_WANTED
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_MULTI',0,ZHOOK_HANDLE)
IPOINTS = SELF%NPOINTS_BLOCK(KBLOCK)
ILEVELS = SELF%NLEVELS(SELF%METADATA(KIDS(1))%D2TYPE)

IF (SELF%MSTORTYPE == JP_STORAGE1D) THEN

  ! Check that the requested fields can be easily assembled into
  ! a 3D array (e.g. without copying): the storage must be contiguous
  ! and the dimensions must be consistent. NFTBC if this kind of access
  ! proves popular and these requirements are not easy to fulfil, we may
  ! need to flexibly "copy in and copy out" if necessary.
  IVARS = 0
  DO I = 1, SIZE(KIDS)
    ID = KIDS(I)
    IF(SELF%FIELD_INDEX(ID)%L_ON) THEN
      IVARS = IVARS+1
      IF (IVARS==1) THEN
        ISTART  = SELF%INDEX_BLOCK(ID,KBLOCK)%ISTART
        IEND    = SELF%INDEX_BLOCK(ID,KBLOCK)%IEND
      ELSE
        IF(SELF%INDEX_BLOCK(ID,KBLOCK)%ISTART == IEND + 1 .AND. &
         & SELF%NLEVELS(SELF%METADATA(ID)%D2TYPE) == ILEVELS) THEN
          IEND = SELF%INDEX_BLOCK(ID,KBLOCK)%IEND
        ELSE
          CALL ABOR1('FIELD_OPEN_MULTI: FIELDS ARE NON-CONTIGUOUS OR DIFFERENT SHAPES: '//SELF%METADATA(ID)%NAME)
          NULLIFY(FPTR)
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  IF (IVARS>= 1) THEN

     FPTR(1:IPOINTS,1:ILEVELS,1:SIZE(KIDS)) => SELF%STORAGE1D(ISTART:IEND)

  ELSE
     NULLIFY(FPTR)
  ENDIF

ELSEIF(SELF%MSTORTYPE == JP_STORAGE_ARR) THEN
  CALL ABOR1('FIELD_CONTAINER_GP_MOD:FMAP_MULTI - NOT IMPLEMENTED FOR JP_STORAGE_ARR')
ELSE

#ifndef FIELD_CONTAINER_GP_MOD_TEST

  ! Unrestricted access into the GFLs - we ignore the list of fids
  CALL GFL_AS_A_WHOLE(KBLOCK, SELF%M_ATTACH, SELF%YRGFL,FPTR)

  ! Special trick for the last, smaller, nproma block: take a copy of the
  ! correct size. Necessary because of the non rank-compliant F77 interfaces
  ! in the IFS (e.g. gprcp etc.). We can't give the user a strided 2D pointer
  IF(KBLOCK == SELF%NBLOCKS .AND. ASSOCIATED(FPTR) .AND. .NOT. SELF%L_ATTACH_GFL_NPROMA_FIXED) THEN

    ISIZE_GFL = SIZE(FPTR,3)
    I2D       = SIZE(SELF%STORAGE_LAST_BLOCK,3)
    IF(I2D < ISIZE_GFL) WRITE(0,*) 'FIELD_OPEN_MULTI: GFL FIELD MAPPING HAS GONE WRONG'
    IF(ILEVELS /= SIZE(FPTR,2)) CALL ABOR1('FIELD_OPEN_MULTI: INCONSISTENT LEVELS')
    SELF%STORAGE_LAST_BLOCK(:,1:ILEVELS,1:MIN(ISIZE_GFL,I2D)) = FPTR(1:IPOINTS,:,1:MIN(ISIZE_GFL,I2D))
    FPTR => SELF%STORAGE_LAST_BLOCK(:,1:ILEVELS,:)

  ENDIF

#endif

ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:FMAP_MULTI',1,ZHOOK_HANDLE)

END SUBROUTINE FMAP_MULTI
!=======================================================================================
SUBROUTINE UNMAP_LAST_BLOCK(SELF, KBLOCK, KID)
! -------------------------------------------------
! GFL-backwards compatibility only: copy last
! block back into GFL. 2D only. Counterpart
! of the "special trick" in fmap_2d
! -------------------------------------------------

CLASS(FIELD_CONTAINER_GP),TARGET, INTENT(IN)  :: SELF
INTEGER(KIND=JPIM),       INTENT(IN)  :: KBLOCK, KID

REAL(KIND=JPRB), POINTER :: FPTR(:,:) => NULL()
INTEGER(KIND=JPIM)       :: IPOINTS, I2D, ILEVELS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:UNMAP_LAST_BLOCK',0,ZHOOK_HANDLE)
IF(SELF%FIELD_INDEX(KID)%L_ON) THEN

  IPOINTS = SELF%NPOINTS_BLOCK(KBLOCK)
  ILEVELS = SELF%NLEVELS(SELF%METADATA(KID)%D2TYPE)

#ifndef FIELD_CONTAINER_GP_MOD_TEST

  CALL GFL_ATTACH(SELF%METADATA,SELF%GFL_MAP_INTERNAL,IPOINTS, ILEVELS, KID, KBLOCK, SELF%M_ATTACH, SELF%YRGMV,SELF%YRGFL, &
                     & SELF%YRSURF,FPTR_2D=FPTR, &
                     & FPTR_UNASSIGNED=SELF%UNASSIGNED_GFL)
  IF(ASSOCIATED(FPTR)) THEN
    IF(.NOT.ASSOCIATED(FPTR, TARGET=SELF%UNASSIGNED_GFL(:,:,1))) THEN

      I2D = SELF%INDEX_LAST_BLOCK(KID)
      FPTR = SELF%STORAGE_LAST_BLOCK(:,1:ILEVELS,I2D)

    ENDIF
  ENDIF

#endif

ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:UNMAP_LAST_BLOCK',1,ZHOOK_HANDLE)

END SUBROUTINE UNMAP_LAST_BLOCK

!=======================================================================================

FUNCTION KINDEX_OLD(SELF,KFID,KATTACH)
! Temporary routine to propagete old GFL index (JGFL)
CLASS(FIELD_CONTAINER_GP) , INTENT(IN) :: SELF
INTEGER(KIND=JPIM),    INTENT(IN) :: KFID
INTEGER(KIND=JPIM),    INTENT(IN) :: KATTACH! Attachment time level (e.g. T9, T0 etc.)
INTEGER(KIND=JPIM) :: KINDEX_OLD
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:KINDEX_OLD',0,ZHOOK_HANDLE)
#ifndef FIELD_CONTAINER_GP_MOD_TEST
KINDEX_OLD=GFL_MAP(SELF%GFL_MAP_INTERNAL,KFID,KATTACH,LD_IGNORE_NOT_SET=.TRUE.)
#endif
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_GP_MOD:KINDEX_OLD',1,ZHOOK_HANDLE)

END FUNCTION KINDEX_OLD
!=======================================================================================

END MODULE FIELD_CONTAINER_GP_MOD
