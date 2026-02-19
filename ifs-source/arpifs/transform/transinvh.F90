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

SUBROUTINE TRANSINVH(YDGEOMETRY,YDGFL,YDGMV,YDML_DIAG,YDML_GCONF,YDDYNA,LDRFRIC,PKRF,CDCONF,YDSP,YDMTRAJ,LDFPOS)

!**** *TRANSINVH * - Inverse transforms

!     Purpose.  Perform inverse transform (spectral to gridpoint)
!     --------

!     Explicit arguments :
!     --------------------
!        CDCONF     - configuration of work
!                     Current values in use for CDCONF(2:3): AA, M1, K0, E0, GB, KD
!        LDFPOS     - post-processing will need the transformed fields

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-10-25
!        R. El Khatib : 01-03-07 Fix
!        P.Smolikova  : 02-09-30 Interface to INI2SCAN2M for d4 in NH
!        R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!        Modified : 03-08-01 M.Hamrud - GFL introduction
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!        K.Yessad (Dec 2003): cleaning in horizontal diffusion.
!        N Wedi (Nov 2011) : add LGRADSP and LLVFP
!        R. El Khatib 17-Jul-2012 Fullpos move away
!        K. Yessad (July 2014): Move some variables.
!        O. Jaron (Sept 2018): introduce LLUVDER for LDFPOS=T for TRANINV_MDL
!     ------------------------------------------------------------------

USE MODEL_DIAGNOSTICS_MOD  , ONLY : MODEL_DIAGNOSTICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE MTRAJ_MOD              , ONLY : MTRAJ
USE YOMGFL                 , ONLY : TGFL
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0                 , ONLY : NCONF, LECMWF
USE YOMSP                  , ONLY : SPVOR_FLT, SPDIV_FLT  
USE YOMSP5                 , ONLY : SPA5
USE YOMDYNA                , ONLY : TDYNA
USE SPECTRAL_FIELDS_MOD
USE YOMLUN                 , ONLY : NULOUT !JARON


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
TYPE(TGMV)                   , INTENT(INOUT) :: YDGMV
TYPE(MODEL_DIAGNOSTICS_TYPE) , INTENT(INOUT) :: YDML_DIAG
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(IN)    :: YDML_GCONF
TYPE(TDYNA)                  , INTENT(IN)    :: YDDYNA
LOGICAL                      , INTENT(IN)    :: LDRFRIC
REAL(KIND=JPRB)              , INTENT(IN)    :: PKRF(:)
CHARACTER(LEN=9)             , INTENT(IN)    :: CDCONF
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSP
TYPE(MTRAJ)                  , INTENT(INOUT), OPTIONAL :: YDMTRAJ
LOGICAL                      , INTENT(IN), OPTIONAL    :: LDFPOS
!     ------------------------------------------------------------------

LOGICAL :: LLDERR,LLUVDER,LLMODE,LLTRAJ,LLFSCOMP,LLVOR,LL_USE_MTRAJ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "transinv_mdl.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSINVH',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDDIMF=>YDML_GCONF%YRDIMF,YDLDDH=>YDML_DIAG%YRLDDH, &
  & YDSPDDH=>YDML_DIAG%YRSPDDH,YDGPDDH=>YDML_DIAG%YRGPDDH)

ASSOCIATE(LADER=>YDDIMF%LADER, LVOR=>YDDIMF%LVOR, LUVDER=>YDDIMF%LUVDER,&
 & GFL=>YDGFL%GFL, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, &
 & GFLTNDHD_DDH=>YDGPDDH%GFLTNDHD_DDH, GMVTNDHD_DDH=>YDGPDDH%GMVTNDHD_DDH, &
 & GMVTNDSI_DDH=>YDGPDDH%GMVTNDSI_DDH, &
 & LRDDHDYN=>YDLDDH%LRDDHDYN, LRHDDDH=>YDLDDH%LRHDDDH, LRSIDDH=>YDLDDH%LRSIDDH, &
 & SPTNDHD_DIV=>YDSPDDH%SPTNDHD_DIV, SPTNDHD_GFL=>YDSPDDH%SPTNDHD_GFL, &
 & SPTNDHD_SNHX=>YDSPDDH%SPTNDHD_SNHX, SPTNDHD_SPD=>YDSPDDH%SPTNDHD_SPD, &
 & SPTNDHD_SVD=>YDSPDDH%SPTNDHD_SVD, SPTNDHD_T=>YDSPDDH%SPTNDHD_T, &
 & SPTNDHD_VOR=>YDSPDDH%SPTNDHD_VOR, SPTNDSI_DIV=>YDSPDDH%SPTNDSI_DIV, &
 & SPTNDSI_SPD=>YDSPDDH%SPTNDSI_SPD, SPTNDSI_SVD=>YDSPDDH%SPTNDSI_SVD, &
 & SPTNDSI_T=>YDSPDDH%SPTNDSI_T, SPTNDSI_VOR=>YDSPDDH%SPTNDSI_VOR)
!     ------------------------------------------------------------------

LL_USE_MTRAJ = PRESENT(YDMTRAJ)

LLDERR = LLE(CDCONF(2:2),'F') .AND. LADER
LLMODE = (CDCONF(2:3)=='AA').OR.(CDCONF(2:3)=='M1').OR.(CDCONF(2:3)=='GB')
LLTRAJ = (CDCONF(2:3)=='E0').OR.(CDCONF(2:3)=='K0').OR.(CDCONF(2:3)=='KD')
LLFSCOMP = LLDERR .AND. (LDRFRIC.AND.YDDYNA%LSLAG.AND.(NCONF == 1.OR.NCONF == 302))
IF (PRESENT(LDFPOS)) THEN
  LLVOR = LDFPOS .OR. (LVOR .AND. LLDERR)
  LLUVDER = (LDFPOS .AND. .NOT. LECMWF) .OR. LUVDER 
  !rek WRITE(UNIT=NULOUT,FMT='(''LDFPOS  = '',L2,'' LLUVDER  = '',L2)') LDFPOS,LLUVDER !TEST
ELSE
  LLVOR = (LVOR .AND. LLDERR)
  LLUVDER = LUVDER
ENDIF

IF(LLMODE) THEN
  IF (YDDYNA%LGRADSP) THEN
    CALL TRANSINV_MDL(YDGEOMETRY,YDGMV,YDML_DIAG%YRMDDH,YDML_GCONF,YDDYNA,YDSP%VOR,YDSP%DIV,YDSP%SP,YDSP%HV,YDSP%GFL,&
     & GMV,GMVS,GFL,LLUVDER,LLDERR,LLVOR,LLFSCOMP,&
     & SPTNDSI_VOR,SPTNDSI_DIV,SPTNDSI_T,&
     & SPTNDSI_SPD,SPTNDSI_SVD,&
     & SPTNDHD_VOR,SPTNDHD_DIV,SPTNDHD_T,&
     & SPTNDHD_SPD,SPTNDHD_SVD,SPTNDHD_SNHX,&
     & SPTNDHD_GFL,&
     & GMVTNDSI_DDH,GMVTNDHD_DDH,GFLTNDHD_DDH,&
     & LRDDHDYN,LRSIDDH,LRHDDDH,PSPVOR_FLT=SPVOR_FLT,PSPDIV_FLT=SPDIV_FLT,PKRF=PKRF)
  ELSE
    CALL TRANSINV_MDL(YDGEOMETRY,YDGMV,YDML_DIAG%YRMDDH,YDML_GCONF,YDDYNA,YDSP%VOR,YDSP%DIV,YDSP%SP,YDSP%HV,YDSP%GFL,&
     & GMV,GMVS,GFL,LLUVDER,LLDERR,LLVOR,LLFSCOMP,&
     & SPTNDSI_VOR,SPTNDSI_DIV,SPTNDSI_T,&
     & SPTNDSI_SPD,SPTNDSI_SVD,&
     & SPTNDHD_VOR,SPTNDHD_DIV,SPTNDHD_T,&
     & SPTNDHD_SPD,SPTNDHD_SVD,SPTNDHD_SNHX,&
     & SPTNDHD_GFL,&
     & GMVTNDSI_DDH,GMVTNDHD_DDH,GFLTNDHD_DDH,&
     & LRDDHDYN,LRSIDDH,LRHDDDH,PKRF=PKRF)
  ENDIF
ELSEIF(LLTRAJ) THEN
  IF (LL_USE_MTRAJ) THEN
    CALL TRANSINV_MDL(YDGEOMETRY,YDGMV,YDML_DIAG%YRMDDH,YDML_GCONF,YDDYNA,SPA5%VOR,SPA5%DIV,SPA5%SP,SPA5%HV,SPA5%GFL,&
     & YDMTRAJ%YRGMV5%GMV5,YDMTRAJ%YRGMV5%GMV5S,YDMTRAJ%YRGFL5%GFL5,LLUVDER,LLDERR,LLVOR,LLFSCOMP,PKRF=PKRF)
  ELSE
    CALL ABOR1('We should not be here YD_TRAJ')
  ENDIF
ELSE
  CALL ABOR1(' TRANSINVH - CONFIGURATION NOT SUPPORTED ')
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TRANSINVH',1,ZHOOK_HANDLE)
END SUBROUTINE TRANSINVH
