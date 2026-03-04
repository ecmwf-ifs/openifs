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

MODULE YOMATLAS

USE PARKIND1, ONLY : JPIM, JPRB
#ifdef WITH_ATLAS
USE ATLAS_MODULE
#endif

IMPLICIT NONE

PRIVATE

PUBLIC :: TATLAS
PUBLIC :: DEALLO_TATLAS

!     ------------------------------------------------------------------

!* ATLAS

TYPE :: TATLAS
#ifdef WITH_ATLAS
  TYPE(ATLAS_GRIDDISTRIBUTION)                :: GRIDDISTRIBUTION
  TYPE(ATLAS_GAUSSIANGRID)                    :: GRID
  TYPE(ATLAS_MESH)                            :: MESH
  TYPE(ATLAS_FVM_METHOD)                      :: FVM
  TYPE(ATLAS_NABLA)                           :: NABLA
  TYPE(ATLAS_FUNCTIONSPACE_NODECOLUMNS)       :: FS_NODECOLUMNS
  TYPE(ATLAS_FUNCTIONSPACE_STRUCTUREDCOLUMNS) :: FS_STRUCTUREDCOLUMNS
#endif
END TYPE TATLAS

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE DEALLO_TATLAS(YDATLAS)

USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

TYPE(TATLAS), POINTER, INTENT(INOUT) :: YDATLAS
#ifdef WITH_ATLAS
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMATLAS:DEALLO_TATLAS',0,ZHOOK_HANDLE)

IF (ASSOCIATED(YDATLAS)) THEN
  CALL YDATLAS%GRID%FINAL()
  CALL YDATLAS%MESH%FINAL()
  CALL YDATLAS%FVM%FINAL()
  CALL YDATLAS%NABLA%FINAL()
  CALL YDATLAS%FS_NODECOLUMNS%FINAL()
ENDIF
NULLIFY(YDATLAS)

IF (LHOOK) CALL DR_HOOK('YOMATLAS:DEALLO_TATLAS',1,ZHOOK_HANDLE)
#endif
END SUBROUTINE DEALLO_TATLAS

! ------------------------------------------------------------------

END MODULE YOMATLAS
