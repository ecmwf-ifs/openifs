! (C) Copyright 1989- Meteo-France.

SUBROUTINE FULLPOS_PRECOND(YDGEOMETRY,YDFIELDS,YDMODEL,YDGMV,YDGFL,KSTEP,KSTOP,YDSURF,PCFUBUF,PXFUBUF)

!**** *FULLPOS_PRECOND*  - Fullpos data preconditioning

!     Purpose.
!     --------
!        To precondition input data prior to Fullpos (ie file packing/unpacking of historical fields)

!**   Interface.
!     ----------
!        *CALL* *FULLPOS_PRECOND(...)

!        Explicit arguments :
!        --------------------
!           YDGEOMETRY : input model geometry
!           YDFIELDS : input model fields
!           KSTEP : current model time step
!           KSTOP : last model time step

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
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPECTRAL_FIELDS_MOD
USE SURFACE_FIELDS_MIX , ONLY : TSURF, CLONE_SURF
USE YOMGMV       , ONLY : TGMV
USE YOMGFL       , ONLY : TGFL

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(GEOMETRY)   ,INTENT(IN) :: YDGEOMETRY
TYPE(FIELDS)      , INTENT(INOUT) :: YDFIELDS
TYPE(MODEL)       , INTENT(IN) :: YDMODEL
TYPE(TGMV)        , INTENT(INOUT) :: YDGMV
TYPE(TGFL)        , INTENT(INOUT) :: YDGFL
INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP
INTEGER(KIND=JPIM), INTENT(IN) :: KSTOP
TYPE(TSURF)       , INTENT(OUT), OPTIONAL :: YDSURF
REAL (KIND=JPRB), ALLOCATABLE, INTENT(OUT), OPTIONAL :: PCFUBUF (:,:,:) ! copy of cumulated fluxes buffer
REAL (KIND=JPRB), ALLOCATABLE, INTENT(OUT), OPTIONAL :: PXFUBUF (:,:,:) ! copy of instantaneous fluxes buffer

!     ------------------------------------------------------------------

LOGICAL :: LLRQCFU, LLRQXFU, LLPACK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "predynfpos.intfb.h"
#include "pregpfpos.intfb.h"

#include "abor1.intfb.h"

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FULLPOS_PRECOND',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,NGPBLKS=>YDDIM%NGPBLKS)
!      -----------------------------------------------------------

IF (.NOT.(PRESENT(YDSURF).AND.PRESENT(PCFUBUF).AND.PRESENT(PXFUBUF))) THEN
  CALL ABOR1('FULLPOS_PRECOND : YDSURF PCFUBUF AND PXFUBUF SHOULD ALL BE PRESENT !')
ENDIF
! Actually, packing should not be needed if these fields are not needed, depending of the post-processing request.
CALL CLONE_SURF(YDGEOMETRY%YRDIM,YDFIELDS%YRSURF,YDSURF)
ALLOCATE(PCFUBUF(NPROMA,YDFIELDS%YRCFU%NFDCFU,NGPBLKS))
ALLOCATE(PXFUBUF(NPROMA,YDFIELDS%YRXFU%NFDXFU,NGPBLKS))
PCFUBUF(:,:,:) = YDFIELDS%YRCFU%GFUBUF(:,:,:)
PXFUBUF(:,:,:) = YDFIELDS%YRXFU%XFUBUF(:,:,:)
LLRQCFU=(YDFIELDS%YRCFU%NFDCFU > 0).AND.(YDFIELDS%YRCFU%LREACFU.OR.KSTEP > 0)
LLRQXFU=(YDFIELDS%YRXFU%NFDXFU > 0).AND.(YDFIELDS%YRXFU%LREAXFU.OR.KSTOP > 0)
CALL PREGPFPOS(YDGEOMETRY,YDSURF,YDFIELDS%YRCFU,YDFIELDS%YRXFU,PXFUBUF,PCFUBUF,LLRQXFU,LLRQCFU)

! GMV/GFL
YDGMV%YT0 = YDFIELDS%YRGMV%YT0
YDGMV%NDIMGMV = YDFIELDS%YRGMV%NDIMGMV
YDGMV%NDIMGMVS = YDFIELDS%YRGMV%NDIMGMVS
ALLOCATE(YDGMV%GMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT0%NDIM,YDGEOMETRY%YRDIM%NGPBLKS))
ALLOCATE(YDGMV%GMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT0%NDIMS,YDGEOMETRY%YRDIM%NGPBLKS))
ALLOCATE(YDGFL%GFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM0,YDGEOMETRY%YRDIM%NGPBLKS))
YDGFL%GFL(:,:,:,:)=YDFIELDS%YRGFL%GFL(:,:,:,:) ! at least to preserve the gp gfl :-(
LLPACK=.TRUE.
CALL PREDYNFPOS(YDGEOMETRY,YDMODEL%YRML_DIAG%YRMDDH,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA,YDFIELDS%YRSPEC,YDGMV,YDGFL,LLPACK)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FULLPOS_PRECOND',1,ZHOOK_HANDLE)
END SUBROUTINE FULLPOS_PRECOND
