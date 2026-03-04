! (C) Copyright 1989- Meteo-France.

SUBROUTINE PREDYNFPOS(YDGEOMETRY,YDMDDH,YDML_GCONF,YDDYNA,YDSPEC,YDGMV,YDGFL,LDPACK)

!**** *PREDYNFPOS*  - Fullpos process

!     Purpose.
!     --------
!        To process Fullpos

!**   Interface.
!     ----------
!        *CALL* *PREDYNFPOS(...)

!        Explicit arguments :
!        --------------------
!           YDGEOMETRY : input model geometry
!           LDPACK     : pack/unpack data before use


!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    
!     ----------    

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib  *METEO-FRANCE*
!        Original : 31-Jul-2012 from PRESPFPOS

!     Modifications.
!     --------------
!        J Hague  : Sep-2013  optional YDSP added
!        K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El khatib 16-May-2014 Optimization of in-line/off-line post-processing reproducibility
!      A. Geer      27-Jul-2015   VarBC is now an object passed by argument, for OOPS
!      M.Hamrud     12-July-2017  Move conversion Tv to T inside LDSPECTRAL, OOPS already T
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYNA                , ONLY : TDYNA
USE YOMMDDH      , ONLY : TMDDH
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV   , ONLY : TGMV
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM
USE SPECTRAL_FIELDS_MOD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN) :: YDGEOMETRY
TYPE(TMDDH)          ,INTENT(IN) :: YDMDDH
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TDYNA)        ,INTENT(IN)    :: YDDYNA
TYPE(SPECTRAL_FIELD) ,INTENT(IN) :: YDSPEC
TYPE(TGMV) ,INTENT(INOUT) :: YDGMV
TYPE(TGFL) ,INTENT(INOUT) :: YDGFL
LOGICAL, INTENT(IN) :: LDPACK

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JKGLO, ICEND, IBL
TYPE(SPECTRAL_FIELD) :: YLSPEC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "transinv_mdl.intfb.h"
#include "etransinv_mdl.intfb.h"
#include "prespfpos.intfb.h"
#include "ctvtot.intfb.h"

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PREDYNFPOS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, NPROMA=>YDDIM%NPROMA)
!      -----------------------------------------------------------

IF (LDPACK) THEN

  ! Work on a copy of the model spectral data
  CALL ALLOCATE_SPEC(YLSPEC,YDSPEC)
  YLSPEC = YDSPEC
  ! Pack/unpack spectral data
  CALL PRESPFPOS(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,YLSPEC)
  ! Inverse transforms
  IF (LELAM) THEN
    CALL ETRANSINV_MDL(YDGEOMETRY,YDGMV,YDMDDH,YDML_GCONF,YLSPEC%VOR,YLSPEC%DIV,&
     & YLSPEC%UB,YLSPEC%VB,&
     & YLSPEC%SP,YLSPEC%HV,YLSPEC%GFL,&
     & YDGMV%GMV,YDGMV%GMVS,YDGFL%GFL,LDERR=.TRUE.)
  ELSE
    CALL TRANSINV_MDL(YDGEOMETRY,YDGMV,YDMDDH,YDML_GCONF,YDDYNA,YLSPEC%VOR,YLSPEC%DIV,&
     & YLSPEC%SP,YLSPEC%HV,YLSPEC%GFL,&
     & YDGMV%GMV,YDGMV%GMVS,YDGFL%GFL,LDUV=.TRUE.,LDERR=.TRUE.,&
     & LDVOR=.TRUE.,LDFSCOMP=.FALSE.)
  ENDIF
  CALL DEALLOCATE_SPEC(YLSPEC)

ELSE

  ! Inverse transforms from the native model data
  IF (LELAM) THEN
    CALL ETRANSINV_MDL(YDGEOMETRY,YDGMV,YDMDDH,YDML_GCONF,YDSPEC%VOR,YDSPEC%DIV,&
     & YDSPEC%UB,YDSPEC%VB,&
     & YDSPEC%SP,YDSPEC%HV,YDSPEC%GFL,&
     & YDGMV%GMV,YDGMV%GMVS,YDGFL%GFL,LDERR=.TRUE.)
  ELSE
    CALL TRANSINV_MDL(YDGEOMETRY,YDGMV,YDMDDH,YDML_GCONF,YDDYNA,YDSPEC%VOR,YDSPEC%DIV,&
     & YDSPEC%SP,YDSPEC%HV,YDSPEC%GFL,&
     & YDGMV%GMV,YDGMV%GMVS,YDGFL%GFL,LDUV=.TRUE.,LDERR=.TRUE.,&
     & LDVOR=.TRUE.,LDFSCOMP=.FALSE.)
  ENDIF

  IF (YDDYNA%LSPRT) THEN
  ! Convert Tv into gridpoint T if (lsprt)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO,ICEND,IBL)
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      CALL CTVTOT(YDGEOMETRY,YDML_GCONF%YGFL,1,ICEND,YDGMV%GMV(1,1,YDGMV%YT0%MT,IBL),YDGFL%GFL(1,1,1,IBL))
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
ENDIF


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PREDYNFPOS',1,ZHOOK_HANDLE)
END SUBROUTINE PREDYNFPOS
