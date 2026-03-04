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

SUBROUTINE TRMFIXERS(YDGEOMETRY,YDGMV,YDCOMPO,YDDPHY,YDML_GCONF,YDML_DYN,LD_CALL_SL,KL0,PLSCAW,PB1,YDSL,PGMVS,PGMVT1S,PGFL,&
                   & PGFLT1,PEXTRADYN)
!
!     Purpose.   INTERFACE SUBROUTINE FOR MASS FIXERS
!     --------   
!
!*   Interface.
!     ----------
!
!        *CALL* *TRMFIXERS(...)*
!
!     Explicit arguments :
!     --------------------
!
!
!!     INPUT:
!     -------------
!        LD_CALL_SL  : true if CALL_SL() subroutine has been called
!        KL0         : indices of the four western points of the 16 point
!                      interpolation grid
!        PLSCAW      : linear weights (distances) for interpolations.
!        PB1         : field to be interpolated
!        YDSL        : SL_STRUCT definition
!        PGMVS       : surface GMV variables at t and t-dt.
!        PGFL        : GFL variables buffer t0 t1
!
!     INPUT/OUTPUT:
!     -------------
!        PGMVT1S      : surface GMV variables at t +dt 
!        PGFLT1       : GFL variables buffer t + dt
!        PEXTRADYN    : mass fixer diagnostics (correction field)
!
!        Implicit arguments :  None.
!        --------------------
!
!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Michail Diamantakis   *ECMWF*

! Modifications
! -------------
! ----------------------------------------------------------------------------

USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMCOMPO               , ONLY : TCOMPO
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE EINT_MOD               , ONLY : SL_STRUCT
USE YOMDPHY                , ONLY : TDPHY

! ----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TCOMPO)      ,INTENT(INOUT) :: YDCOMPO
TYPE(TDPHY)       ,INTENT(INOUT) :: YDDPHY
TYPE(MODEL_DYNAMICS_TYPE)    ,INTENT(INOUT):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
LOGICAL           ,INTENT(IN)    :: LD_CALL_SL
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDGEOMETRY%YRDIM%NGPBLKS)
TYPE(SL_STRUCT)   ,INTENT(IN)    :: YDSL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM, &
 & YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM1, &
 & YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEXTRADYN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDDPHY%NVEXTRDYN,&
 & YDGEOMETRY%YRDIM%NGPBLKS)
! ----------------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "negfixer.intfb.h"
#include "qmfixer.intfb.h"
#include "qmfixer2.intfb.h"
#include "jmgfixer.intfb.h"
#include "tracmf.intfb.h"

IF (LHOOK) CALL DR_HOOK('TRMFIXERS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(LADVNEGFIX=>YGFL%LADVNEGFIX, LTRCMFBC=>YGFL%LTRCMFBC, &
 & LTRCMFIX=>YGFL%LTRCMFIX, LTRCMFMG=>YGFL%LTRCMFMG, LTRCMFP=>YGFL%LTRCMFP, &
 & LTRCMFPR=>YGFL%LTRCMFPR, NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, &
 & YCOMP=>YGFL%YCOMP, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NVEXTRDYN=>YDDPHY%NVEXTRDYN, &
 & NDIMGMVS=>YDGMV%NDIMGMVS, YT1=>YDGMV%YT1, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1)
! ----------------------------------------------------------------------------

! Apply simple negative fixer if required 
IF (LADVNEGFIX) CALL NEGFIXER(YDGEOMETRY,YDGMV,YDML_GCONF,YDPTRSLB1,KL0,PB1,YDSL,PGMVS,PGMVT1S,PGFL,PGFLT1)

IF (LTRCMFIX) THEN
  IF (LTRCMFBC.OR.LTRCMFP) THEN
    ! Bermejo-Conde/ proportional mass fixer
    IF (.NOT.LD_CALL_SL) THEN
      CALL ABOR1(' TRMFIXERS: QMFIXER requires output of CALL_SL!')
    ENDIF
    CALL QMFIXER(YDGEOMETRY,YDGMV,YDCOMPO,YDDPHY,YGFL,YDML_DYN,KL0,PLSCAW,PB1,YDSL,PGMVS,PGMVT1S,PGFL,PGFLT1,PEXTRADYN)
  ELSEIF (LTRCMFPR) THEN
    ! Priestley mass fixer
    IF (.NOT.LD_CALL_SL) THEN
      CALL ABOR1(' TRMFIXERS: QMFIXER2 requires output of CALL_SL!')
    ENDIF
    CALL QMFIXER2(YDGEOMETRY,YDGMV,YDDPHY,YGFL,YDML_DYN,KL0,PLSCAW,PB1,YDSL,PGMVS,PGMVT1S,PGFL,PGFLT1,PEXTRADYN)
  ELSEIF (LTRCMFMG) THEN
    ! Mc Gregor's mass fixer
    IF (.NOT.LD_CALL_SL) THEN
      CALL ABOR1(' TRMFIXERS: JMGFIXER requires output of CALL_SL!')
    ENDIF
    CALL JMGFIXER(YDGEOMETRY,YDGMV,YDDPHY,YGFL,YDML_DYN%YRDYN,YDPTRSLB1,KL0,PB1,YDSL,PGMVS,PGMVT1S,PGFL,PGFLT1,PEXTRADYN)
  ELSE
    ! additive mass fixers
    CALL TRACMF(YDGEOMETRY,YDGMV,YGFL,PGMVS,PGMVT1S,PGFL,PGFLT1)
  ENDIF
ENDIF

!* END PART FOR  TRMFIXERS

! ----------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TRMFIXERS',1,ZHOOK_HANDLE)
END SUBROUTINE TRMFIXERS
