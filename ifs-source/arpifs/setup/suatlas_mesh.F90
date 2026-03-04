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

SUBROUTINE SUATLAS_MESH(YDGEOMETRY)

!**** *SUATLAS_MESH*  - SETUP ATLAS MESH ETC.

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        *CALL* *SUATLAS_MESH*

!        EXPLICIT ARGUMENTS
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Mats Hamrud AND Willem Deconinck  *ECMWF*
!      ORIGINAL : Nov 2015

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
#ifdef WITH_ATLAS
USE PARKIND1 , ONLY : JPIM,JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMCST   , ONLY : RA
USE YOM_ATLAS_IFS, ONLY: LATLAS_MESH
USE ATLAS_MODULE, ONLY: ATLAS_GRIDDISTRIBUTION, ATLAS_MESHGENERATOR, ATLAS_FUNCTIONSPACE_STRUCTUREDCOLUMNS, ATLAS_NABLA, &
  & ATLAS_CONFIG, ATLAS_REDUCEDGAUSSIANGRID, ATLAS_FVM_METHOD
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
#ifdef WITH_ATLAS

! Objects that are temporary during initialisation
TYPE(ATLAS_GRIDDISTRIBUTION) :: DISTRIBUTION
TYPE(ATLAS_MESHGENERATOR)    :: MESHGENERATOR
TYPE(ATLAS_CONFIG)           :: CONFIG

INTEGER(KIND=JPIM) :: IFVM_HALO_SIZE = 2
INTEGER(KIND=JPIM) :: ISTRUCTURED_HALO_SIZE = 3

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

NAMELIST/NAM_ATLAS_IFS/LATLAS_MESH

#include "posnam.intfb.h"

!     ------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('SUATLAS_MESH',0,ZHOOK_HANDLE)
CALL POSNAM(NULNAM,'NAM_ATLAS_IFS')
READ(NULNAM,NAM_ATLAS_IFS)

IF(LATLAS_MESH)THEN

  WRITE(NULOUT,*) 'ATLAS IS SWITCHED ON. LATLAS_MESH FLAG IS:',LATLAS_MESH
  IF( .NOT. ASSOCIATED(YDGEOMETRY%YRATLAS) ) THEN
    ALLOCATE( YDGEOMETRY%YRATLAS )
  ENDIF
  ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
   & YDATLAS=>YDGEOMETRY%YRATLAS)
  ASSOCIATE(NLOENG=>YDGEM%NLOENG, NGLOBALPROC=>YDMP%NGLOBALPROC, NDGLG=>YDDIM%NDGLG)

  !     ------------------------------------------------------------------

  ! Create a new Reduced Gaussian Grid based on the Gaussian N and the nloen array
  !YDATLAS%GRID = ATLAS_REDUCEDGAUSSIANGRID( NDGLG/2, NLOENG(1:NDGLG/2) )
  YDATLAS%GRID = ATLAS_REDUCEDGAUSSIANGRID( NLOENG(1:NDGLG) )
  IF( YDATLAS%GRID%SIZE() /= YDGEOMETRY%YRGEM%NGPTOTG ) THEN
    CALL ABOR1("Mismatch in Atlas GRID%SIZE() and NGPTOTG")
  ENDIF

  ! Grid distribution: Impose grid distribution from IFS, and partition numbering starts at 1
  YDATLAS%GRIDDISTRIBUTION = ATLAS_GRIDDISTRIBUTION( NGLOBALPROC, PART0=1 )

  ! Generate mesh with given grid and distribution
  CONFIG = ATLAS_CONFIG()
  call CONFIG%SET('triangulate',.FALSE.)
  call CONFIG%SET('angle',0)
  MESHGENERATOR = ATLAS_MESHGENERATOR(CONFIG)
  YDATLAS%MESH  = MESHGENERATOR%GENERATE(YDATLAS%GRID,YDATLAS%GRIDDISTRIBUTION)
  CALL DISTRIBUTION%FINAL()
  CALL MESHGENERATOR%FINAL()
  ! Generate function-space, with a given halo_size
  CONFIG = ATLAS_CONFIG()
  CALL CONFIG%SET('halo',IFVM_HALO_SIZE)
  CALL CONFIG%SET('radius',RA)
  YDATLAS%FVM = ATLAS_FVM_METHOD(YDATLAS%MESH,CONFIG)
  CALL CONFIG%FINAL()
  ! Setup nabla operator
  YDATLAS%NABLA = ATLAS_NABLA(YDATLAS%FVM)
  YDATLAS%FS_NODECOLUMNS = YDATLAS%FVM%NODE_COLUMNS()
  YDATLAS%FS_STRUCTUREDCOLUMNS = ATLAS_FUNCTIONSPACE_STRUCTUREDCOLUMNS( &
    & YDATLAS%GRID, YDATLAS%GRIDDISTRIBUTION, HALO=ISTRUCTURED_HALO_SIZE )

  !     ------------------------------------------------------------------

  END ASSOCIATE
  END ASSOCIATE
ENDIF
IF (LHOOK) CALL DR_HOOK('SUATLAS_MESH',1,ZHOOK_HANDLE)
#endif
END SUBROUTINE SUATLAS_MESH
