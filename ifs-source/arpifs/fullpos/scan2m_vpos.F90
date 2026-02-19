! (C) Copyright 1989- Meteo-France.

SUBROUTINE SCAN2M_VPOS(YDQTYPE,YDNAMFPSCI,YDAFN,LDHPOS,YDVAB,YDGEOMETRY,YDGMV,YDGFL,YDSURF,YDXFU,KCUFNR,PMCUFGP, &
 & YDML_GCONF,YDDYN,YDDYNA,YDML_PHY_MF,CDCONF,KPXLEV,PXLEV,YDTFP_DYNDS,PAUX,PGPP)

!****-------------------------------------------------------------------
!**** *SCAN2M_VPOS* - Fullpos grid-point space computations
!****-------------------------------------------------------------------
!     Purpose.   post-processing in grid-point space
!     --------   

!**   Interface.
!     ----------
!        *CALL* *SCAN2M_VPOS (..)

!        Explicit arguments :  CDCONF - configuration of work (see doc.)
!        --------------------  
!        YDTFP_DYNDS : overall fields descriptors

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Ryad El Khatib *Meteo-France* 

! Modifications
! -------------
!   original 18-Jul-2012 from SCAN2M
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      R. El Khatib & D. Salmond Remove extra GSTATS calls.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!-----------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE YOMXFU       , ONLY : TXFU
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE YOMGFL       , ONLY : TGFL
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LELAM
USE TYPE_GFLFLDS , ONLY : TYPE_IGFLFLDD, TYPE_LGFLFLDD
USE TYPE_GMVS    , ONLY : TYPE_T0
USE YOMVERT      , ONLY : TVAB
USE FULLPOS_MIX  , ONLY : FULLPOS_TYPE
USE YOMAFN       , ONLY : TAFN
USE YOMFPC       , ONLY : TNAMFPSCI
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TYPE_FPRQDYN), INTENT(IN) :: YDQTYPE
TYPE(TNAMFPSCI),    INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN),         INTENT(IN) :: YDAFN
LOGICAL,            INTENT(IN) :: LDHPOS
TYPE(TVAB),         INTENT(IN) :: YDVAB
TYPE(GEOMETRY),     INTENT(IN) :: YDGEOMETRY
TYPE(TGMV),         INTENT(IN) :: YDGMV
TYPE(TGFL),         INTENT(IN) :: YDGFL
TYPE(TSURF),        INTENT(IN) :: YDSURF
TYPE(TDYN),         INTENT(IN) :: YDDYN
TYPE(TDYNA),        INTENT(IN) :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
TYPE(TXFU),         INTENT(IN) :: YDXFU
INTEGER(KIND=JPIM), INTENT(IN) :: KCUFNR
REAL(KIND=JPRB),    INTENT(IN) :: PMCUFGP(YDGEOMETRY%YRDIM%NPROMA,KCUFNR,YDGEOMETRY%YRDIM%NGPBLKS)
CHARACTER(LEN=1),   INTENT(IN) :: CDCONF
INTEGER(KIND=JPIM), INTENT(IN) :: KPXLEV
REAL(KIND=JPRB),    INTENT(IN) :: PXLEV(KPXLEV)
TYPE(FULLPOS_TYPE), INTENT(IN) :: YDTFP_DYNDS(:)
REAL(KIND=JPRB),    INTENT(OUT) :: PAUX(YDGEOMETRY%YRDIM%NPROMA,YDQTYPE%NFPAUXB,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB),    INTENT(OUT) :: PGPP(YDGEOMETRY%YRDIM%NPROMA,YDQTYPE%NFPGT1,YDGEOMETRY%YRDIM%NGPBLKS)

!     ------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IBL, IEND1, IOFF, IST1, JKGLO

TYPE(TYPE_IGFLFLDD) :: YLIN_GFL
TYPE(TYPE_LGFLFLDD) :: YLGFL
TYPE(TYPE_T0)       :: YLIN_GMV

!     ------------------------------------------------------------------


#include "ebipos.intfb.h"
#include "vpos.intfb.h"
#include "vpos_prep.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SCAN2M_VPOS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NDIM=>YGFL%NDIM, &
 & NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, NDIMGMV=>YDGMV%NDIMGMV, &
 & NDIMGMVS=>YDGMV%NDIMGMVS)

!*       2.    INITIAL GRIDPOINT COMPUTATIONS
!              ------------------------------

CALL VPOS_PREP(YDGMV,YGFL,YLIN_GFL,YLGFL,YLIN_GMV)

CALL GSTATS(1435,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IBL,IOFF,IST1,IEND1)
DO JKGLO=1,NGPTOT,NPROMA
  IBL=(JKGLO-1)/NPROMA+1
  IOFF=JKGLO
  IST1=1
  IEND1=MIN(NPROMA,NGPTOT-JKGLO+1)
  CALL VPOS(YDQTYPE,YDNAMFPSCI,YDAFN,LDHPOS,YDVAB,YDGEOMETRY,YDSURF,YDXFU,KCUFNR,PMCUFGP(:,:,IBL),YDML_GCONF,YDDYN,YDDYNA,YDML_PHY_MF, &
   & IST1,IEND1,IOFF,NDIM,NDIMGMV,NDIMGMVS,CDCONF,KPXLEV,PXLEV,YLIN_GFL,YLGFL, &
   & YLIN_GMV, GMV(:,:,:,IBL),GMVS(:,:,IBL),GFL(:,:,:,IBL),YDQTYPE%NFPGT1,YDQTYPE%NFPAUXB,PAUX(:,:,IBL),PGPP(:,:,IBL))
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1435,1)
IF (LELAM) THEN
  CALL EBIPOS(YDQTYPE,YDGEOMETRY,YDQTYPE%NFPGT1,YDTFP_DYNDS(:)%LLBIP,PGPP)
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SCAN2M_VPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SCAN2M_VPOS
