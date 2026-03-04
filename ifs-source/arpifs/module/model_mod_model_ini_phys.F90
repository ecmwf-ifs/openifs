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

SUBMODULE (MODEL_MOD) MODEL_MOD_MODEL_INI_PHYS
IMPLICIT NONE

CONTAINS

MODULE SUBROUTINE MODEL_INI_PHYS(SELF, YDFIELDS, YDTIME)


!**** *MODEL_INI_PHYS*  - Initialize physical diagnostics (OOPS interface to STEPX)

!     Purpose.
!     --------
!        Initialize physical diagnostics

!**   Interface.
!     ----------
!        *CALL* *MODEL_INI_PHYS*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    see includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Ryad El Khatib *Meteo-FRance*
!      Original : 09-Jul-2018 from CNT4

!     Modifications
!     -------------
!      6-Aug-2021  M.Chrust Revert surface fields as STEPX calls grid point computations
!     20-Aug-2021  M.Chrust Call STEPO_OOPS in physics initialization mode rather than STEPX

! Modifications
! -------------
! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1        , ONLY : JPRB, JPIM
USE YOMHOOK         , ONLY : LHOOK, DR_HOOK, JPHOOK
USE TYPE_MODEL      , ONLY : MODEL
USE GEOMETRY_MOD    , ONLY : GEOMETRY_SAME
USE FIELDS_MOD      , ONLY : FIELDS
USE YOMCT0          , ONLY : LCALLSFX
USE YOMCT3          , ONLY : NSTEP
USE YOMGFL          , ONLY : COPY_YOMGFL
USE YOMGMV          , ONLY : COPY_YOMGMV
USE YOMSIG          , ONLY : RESTART
USE YOMTRAJ_OOPS    , ONLY : TRAJ_TYPE_OOPS
USE MTRAJ_MOD       , ONLY : MTRAJ
USE DATETIME_TMP_MOD, ONLY : DATETIME_TMP, SETRIP0
USE SURFACE_FIELDS_MIX, ONLY : TSURF, COPY_SURF, COPY_CTOR_SURF

IMPLICIT NONE

TYPE(MODEL)       , INTENT(INOUT) :: SELF
TYPE(FIELDS)      , INTENT(INOUT) :: YDFIELDS
TYPE(DATETIME_TMP), INTENT(IN)    :: YDTIME

TYPE(MTRAJ)          :: YLMTRAJ
TYPE(TRAJ_TYPE_OOPS) :: YLOOPSTRAJ
TYPE(TSURF)          :: YLSURF

CHARACTER(LEN=9)   :: CL='0AAX00000'
INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPRB), ALLOCATABLE :: ZGFL(:,:,:,:), ZGMV(:,:,:,:), ZGMVS(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "ecradfr.intfb.h"
#include "updtim.intfb.h"




!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MODEL_MOD:MODEL_INI_PHYS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDFIELDS%GEOM%YRDIM, YDDIMV=>YDFIELDS%GEOM%YRDIMV,&
 & YGFL=>YDFIELDS%STATE_MODEL%YRML_GCONF%YGFL, YDGMV=>YDFIELDS%YRGMV,&
 & YDRIP=>SELF%YRML_GCONF%YRRIP,YDDYNA=>SELF%YRML_DYN%YRDYNA)
ASSOCIATE(NDIM=>YGFL%NDIM, NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA,&
 & NFLEVG=>YDDIMV%NFLEVG, NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS)

CALL MODEL_SET(SELF)

IF (.NOT.GEOMETRY_SAME(SELF%YRML_GCONF%GEOM,YDFIELDS%GEOM)) THEN
  CALL ABOR1('MODEL_MOD:MODEL_INI_PHYS - GEOMETRIES DIFFERENT')
ENDIF

CALL SETRIP0(YDTIME)
CALL UPDTIM(YDFIELDS%GEOM,YDFIELDS%YRSURF,SELF,NSTEP, &
 &          SELF%YRML_GCONF%YRRIP%TDT,SELF%YRML_GCONF%YRRIP%TSTEP,.TRUE.)
IF (SELF%YRML_PHY_EC%YREPHY%LEPHYS.OR. &
 & (SELF%YRML_PHY_MF%YRPHY%LMPHYS.AND.SELF%YRML_PHY_MF%YRPHY%LRAYFM)) THEN
  CALL ECRADFR(YDDIM,SELF%YRML_PHY_RAD,YDDYNA,YDRIP)
ENDIF
IF(SELF%YRML_PHY_MF%YRPHY%LMPHYS.AND.SELF%YRML_PHY_MF%YRPHY%LRAYFM15) THEN
  CALL ABOR1('Radiation code from cycle 15 no longer available.')
 !!CALL ECRADFR15(YDDIM,YDRIP)
ENDIF
RESTART='0'
YDFIELDS%YRXFU%XFUBUF(:,:,:)=0.0_JPRB ! at least to fill the extension zone (LAM models)

! Work on copies of GMV/GFL/SURF fields
! because they will be modified by gridpoint computations
ALLOCATE(ZGMV(NPROMA,NFLEVG,NDIMGMV,NGPBLKS))
ALLOCATE(ZGMVS(NPROMA,NDIMGMVS,NGPBLKS))
ALLOCATE(ZGFL(NPROMA,NFLEVG,NDIM,NGPBLKS))
!$OMP PARALLEL DO PRIVATE(JBL)
  DO JBL = 1,NGPBLKS
    ZGMV(:,:,:,JBL) = YDFIELDS%YRGMV%GMV(:,:,:,JBL)
    ZGMVS (:,:,JBL) = YDFIELDS%YRGMV%GMVS (:,:,JBL)
    ZGFL(:,:,:,JBL) = YDFIELDS%YRGFL%GFL(:,:,:,JBL)
  ENDDO
!$OMP END PARALLEL DO
CALL COPY_CTOR_SURF(YLSURF,YDFIELDS%YRSURF)





! Make sure Physics will be called at least once for NSTEP=0 (SURFEX)
LCALLSFX=.TRUE.  ! Return to normal step configuration

! Recover GMV/GFL/SURF fields
!$OMP PARALLEL DO PRIVATE(JBL)
  DO JBL = 1,NGPBLKS
    YDFIELDS%YRGMV%GMV(:,:,:,JBL) = ZGMV(:,:,:,JBL)
    YDFIELDS%YRGMV%GMVS (:,:,JBL) = ZGMVS (:,:,JBL)
    YDFIELDS%YRGFL%GFL(:,:,:,JBL) = ZGFL(:,:,:,JBL)
  ENDDO
!$OMP END PARALLEL DO
CALL COPY_SURF(YDFIELDS%YRSURF,YLSURF)

IF (YDFIELDS%YRCFU%LCUMFU) THEN
  YDFIELDS%YRCFU%GFUBUF(:,:,:)=0.0_JPRB
ENDIF

CALL MODEL_UNSET(SELF)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MODEL_MOD:MODEL_INI_PHYS',1,ZHOOK_HANDLE)

END SUBROUTINE MODEL_INI_PHYS

END SUBMODULE
