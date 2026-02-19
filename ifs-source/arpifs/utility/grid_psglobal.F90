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

SUBROUTINE GRID_PSGLOBAL(YDGEOMETRY,PSPSP)

!     Purpose.
!     --------
!       Compute global ps grid field from distributed lnps spectral field.
!       High resolution gridpoint surface pressure is required for 
!       conserving interpolation.

!     Arguments.
!     ----------
!       PSPSP   : Local lnPs spectral field

!     Author.
!     -------
!       E.Holm (part of PRE_GRID_BICONSERV by Y. Tremolet)

!     Modifications.
!     --------------
!       Original : 12-Feb-2008
! ----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1,ONLY: JPIM, JPRB
USE YOMHOOK ,ONLY: LHOOK, DR_HOOK, JPHOOK
USE YOM_GRID_BICONSERV, ONLY: RGPPRS_HR

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB), INTENT(IN) :: PSPSP(YDGEOMETRY%YRDIM%NSPEC2)

REAL(KIND=JPRB), ALLOCATABLE :: ZPRSGLO(:,:),ZPRSLOC(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gathergpf.intfb.h"
#include "speree.intfb.h"
! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRID_PSGLOBAL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG)

ALLOCATE(ZPRSGLO(NGPTOTG,1))
ALLOCATE(ZPRSLOC(NGPTOT))

!-----------------------------------------------------------------------
! 1. Transform HR lnPs to local gridpoint space and convert from 
!    ln(ps) to ps.
!-----------------------------------------------------------------------
CALL SPEREE(YDGEOMETRY,1,1,PSPSP,ZPRSLOC)
ZPRSLOC(:)=EXP(ZPRSLOC(:))

!-----------------------------------------------------------------------
! 2. Gather local grid point sections to a global field and 
!    write to module
!-----------------------------------------------------------------------
CALL GATHERGPF(YDGEOMETRY,ZPRSLOC,ZPRSGLO,1,-1)

RGPPRS_HR(:)=ZPRSGLO(:,1)

DEALLOCATE(ZPRSGLO)
DEALLOCATE(ZPRSLOC)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GRID_PSGLOBAL',1,ZHOOK_HANDLE)
END SUBROUTINE GRID_PSGLOBAL
