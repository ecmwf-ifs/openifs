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

SUBROUTINE STEPX(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL)

!**** *STEPX*  - Step to initialize physical diagnostics

!     Purpose.
!     --------
!        Initialize physical diagnostics

!**   Interface.
!     ----------
!        *CALL* *STEPX()

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
!      6-Aug-2021  M.Chrust Call grid point computations for OOPS (required to fill in MWAVE fields)
!     20-Aug-2021  M.Chrust Delagated to STEPO_OOPS with LINIPHYSONLY set to true

! Modifications
! -------------
! End Modifications
!------------------------------------------------------------------------------

USE TYPE_MODEL          , ONLY : MODEL
USE GEOMETRY_MOD        , ONLY : GEOMETRY
USE FIELDS_MOD          , ONLY : FIELDS
USE MTRAJ_MOD           , ONLY : MTRAJ
USE PARKIND1            , ONLY : JPRD, JPIM, JPRB
USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN              , ONLY : NULOUT
USE YOMTIM              , ONLY : RSTART, RVSTART, RTIMEF
USE YOMCT0              , ONLY : LCALLSFX, LOPDIS
USE YOMCT3              , ONLY : NSTEP
USE YOMTRAJ_OOPS        , ONLY : TRAJ_TYPE_OOPS
USE YOMSIG              , ONLY : RESTART

!      -----------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT(IN) :: YDGEOMETRY
TYPE(FIELDS)        , INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)         , INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)         , INTENT(INOUT) :: YDMODEL

CHARACTER(LEN=9)   :: CLCONF
INTEGER(KIND=JPIM) :: IJUM
TYPE(TRAJ_TYPE_OOPS) :: YLOOPSTRAJ

LOGICAL :: LL_TST_GPGFL ! do timestepping on grid-point YDFIELDS%YRGFL%GFL
LOGICAL :: LL_DFISTEP

REAL(KIND=JPRD) :: ZCT, ZVT, ZWT, ZJUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------

#include "user_clock.intfb.h"

#include "opdis.intfb.h"
#include "updtim.intfb.h"
#include "updnemoocean.intfb.h"
#include "ecradfr.intfb.h"
#include "zeroddh.intfb.h"




!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('STEPX',0,ZHOOK_HANDLE)
ASSOCIATE(YDGMV5=>YDMTRAJ%YRGMV5, YDSPEC=>YDFIELDS%YRSPEC, &
 & YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NGPTOT=>YDGEM%NGPTOT, NSTART=>YDRIP%NSTART, NSTOP=>YDRIP%NSTOP, &
 & LHDOUFD=>YDLDDH%LHDOUFD, LHDOUFG=>YDLDDH%LHDOUFG, LHDOUFZ=>YDLDDH%LHDOUFZ, LHDOUP=>YDLDDH%LHDOUP)
!      -----------------------------------------------------------

IF (YDFIELDS%YRXFU%LXFU) THEN

!*       1.    VARIOUS INITIALIZATION
!              ----------------------

  CLCONF(1:9)='0AAX00000'

  CALL USER_CLOCK(PELAPSED_TIME=ZWT,PTOTAL_CP=ZCT,PVECTOR_CP=ZVT)
  ZCT=ZCT-RSTART
  ZVT=ZVT-RVSTART
  ZWT=ZWT-RTIMEF
  RSTART=RSTART+ZCT
  RVSTART=RVSTART+ZVT
  RTIMEF=RTIMEF+ZWT
  ZJUM=10._JPRB**(INT(LOG10(REAL(MAX(NSTOP-NSTART,1),JPRB)))+1)
  IJUM=NINT(MAX(ZJUM/100._JPRB,1.0_JPRB))
  IF(NSTEP-NSTART <= 10.OR.MOD(NSTEP,IJUM) == 0)THEN
    WRITE(NULOUT,'('' NSTEP ='',I6,'' STEPX   '',A9)') NSTEP,CLCONF
  ENDIF
  IF (LOPDIS) THEN
    CALL OPDIS(CLCONF,'STEPX   ',ZCT,ZVT,ZWT,RSTART,RTIMEF,NSTEP,YDDYNA%LNHDYN)
  ENDIF

  CALL UPDTIM(YDGEOMETRY,YDFIELDS%YRSURF,YDMODEL,0,YDRIP%TSTEP,YDRIP%TSTEP,.FALSE.)
#ifdef WITH_NEMO
  IF(YDMODEL%YRML_AOC%YRMCC%LMCC04.AND.YDMODEL%YRML_AOC%YRMCC%LNEMOCOUP) THEN
    CALL GSTATS(33,0)
    CALL UPDNEMOOCEAN(YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_EC%YREPHY,YDRIP,NPROMA,NGPTOT,0,YDRIP%TSTEP,YDFIELDS%YRSURF,YDGEOMETRY)
    CALL GSTATS(33,1)
  ENDIF
#endif
  IF (YDMODEL%YRML_PHY_EC%YREPHY%LEPHYS.OR. &
   & (YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS.AND.YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM)) THEN
    CALL ECRADFR(YDDIM,YDMODEL%YRML_PHY_RAD,YDDYNA,YDRIP)
  ENDIF
  IF(YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS.AND.YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM15) THEN
    CALL ABOR1('Radiation code from cycle 15 no longer available.')
    !!CALL ECRADFR15(YDDIM,YDRIP)
  ENDIF
  RESTART='0'
  YDFIELDS%YRXFU%XFUBUF(:,:,:)=0.0_JPRB ! at least to fill the extension zone (LAM models)


!*       2.   GRIDPOINT COMPUTATIONS.
!             -----------------------

  IF (LHDOUFG.OR.LHDOUFZ.OR.LHDOUFD.OR.LHDOUP) THEN
    CALL ZERODDH(YDDIMV,YDMODEL%YRML_DIAG)
  ENDIF

  ! Timestepping on grid-point GFL is done when both:
  ! - grid-point calculations are required
  ! - direct spectral transforms are required
  LL_TST_GPGFL=.FALSE.
  LL_DFISTEP=.FALSE.



  ! Make sure Physics will be called at least once for NSTEP=0 (SURFEX)
  LCALLSFX=.TRUE.  ! Return to normal step configuration

ENDIF

IF (YDFIELDS%YRCFU%LCUMFU) THEN

  YDFIELDS%YRCFU%GFUBUF(:,:,:)=0.0_JPRB

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('STEPX',1,ZHOOK_HANDLE)
END SUBROUTINE STEPX
