! (C) Copyright 1989- Meteo-France.

#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE VPOS(YDQTYPE,YDNAMFPSCI,YDAFN,LDHPOS,YDFPVAB,YDGEOMETRY,YDSURF,YDXFU,KCUFNR,PMCUF,YDML_GCONF,YDDYN,YDDYNA,YDML_PHY_MF, &
 & KST,KEND,KSTGLO,KGFL,KGMV,KGMVS,CDCONF,KPXLEV,PXLEV,YDIN_GFL,YDGFL,YDIN_GMV,PGMV,PGMVS,PGFL,KGT1,KAUX,PAUX,PGPP)

!**** *VPOS*  - VERTICAL POST-PROCESSING - FULL POS

! PURPOSE.
! --------
!  Full-POS interface to POS.

! INTERFACE.
! ----------
!   *CALL* *VPOS*

! EXPLICIT ARGUMENTS
! --------------------
! LDHPOS       : .TRUE. if any horizontal gridpoint post-processing (core or E-zone)
! KST,KEND     : First and last adress of computation in NPROMA-sized arrays
! KSTGLO       : Start adress of computation in NGPTOT-sized arrays
! KGFL         : Number of GFL fields in PGFL
! KGMV         : Number of GMV fields in PGMV
! KGMVS        : Number of GMVS fields in PGMVS
! CDCONF       : Configuration of work (used in POS to know the vertical coord of interpolation)
! YDIN_GFL     : Pointers allowing to retrieve individual fields in PGFL
! YDGFL        : says if GFL is active or not
! YDIN_GMV     : Pointers allowing to retrieve individual fields in PGMV and PGMVS
! PGMV         : GMV fields
! PGMVS        : GMVS fields
! PGFL         : GFL fields
! PAUX         : output data to remain not spectrally fitted
! PGPP         : output data to be spectrally fitted
! PMCUF        : monitoring of coupling update frequency 

! IMPLICIT ARGUMENTS
! --------------------

! METHOD.
! -------
!   SEE DOCUMENTATION

! EXTERNALS.
! ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad (Aug 2009): rewrite vertical interpolator in post-processing.
!      H. Hersbach  : 01-04-2011 auxiliary diagnostic radiation fields
!      K. Yessad (Dec 2011): various contributions.
!      R. El Khatib : 23-Mar-2012 Fix bounds checking issues
!      R. El Khatib 20-Aug-2012 GAUXBUF removed and replaced by HFPBUF
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      K. Yessad (June 2017): Introduce NHQE model.
!      R. El Khatib 11-Oct-2017 Optimization (use pointers rather than copies)
!      E.Dutra/G.Arduini Jan 2018: change of SP_SG to 4 Dimensions, snow multi-layer 
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!      R. El Khatib 14-May-2018 fix a memory leak
!      09-2018 R. Brozkova, A. Bucanek: MOCON diagnostics in offline fullpos
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD   , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOMDYN                 , ONLY : TDYN
USE YOMDYNA                , ONLY : TDYNA
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX     , ONLY : TSURF
USE PARKIND1               , ONLY : JPIM     ,JPRB
USE YOMHOOK                , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMCST                 , ONLY : RG
USE YOMSTA                 , ONLY : RDTDZ1,NLEXTRAP
USE YOMXFU                 , ONLY : TXFU
USE TYPE_GFLFLDS           , ONLY : TYPE_IGFLFLDD,TYPE_LGFLFLDD
USE TYPE_GMVS              , ONLY : TYPE_T0
USE INTDYN_MOD             , ONLY : YYTXYB
USE YOMVERT                , ONLY : TVAB
USE PARDIM                 , ONLY : JPNPPM
USE YOMAFN                 , ONLY : TAFN
USE YOMFPC                 , ONLY : TNAMFPSCI
USE TYPE_FPRQDYNS          , ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(TNAMFPSCI)    ,INTENT(IN)    :: YDNAMFPSCI
TYPE(TAFN)         ,INTENT(IN)    :: YDAFN
LOGICAL            ,INTENT(IN)    :: LDHPOS
TYPE(TVAB)         ,INTENT(IN)    :: YDFPVAB
TYPE(GEOMETRY)     ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)        ,INTENT(IN)    :: YDSURF
TYPE(TDYN)         ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)        ,INTENT(IN)    :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN)  :: YDML_GCONF
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN)  :: YDML_PHY_MF
TYPE(TXFU)         ,INTENT(IN), TARGET    :: YDXFU
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KCUFNR
REAL(KIND=JPRB)    ,INTENT(IN)    :: PMCUF(YDGEOMETRY%YRDIM%NPROMA,KCUFNR)
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KST 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KSTGLO 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KGFL
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KGMV
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KGMVS
CHARACTER(LEN=1)   ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KPXLEV
REAL(KIND=JPRB)    ,INTENT(IN)    :: PXLEV(KPXLEV)
TYPE(TYPE_IGFLFLDD),INTENT(IN)    :: YDIN_GFL
TYPE(TYPE_LGFLFLDD),INTENT(IN)    :: YDGFL
TYPE(TYPE_T0)      ,INTENT(IN)    :: YDIN_GMV
REAL(KIND=JPRB)    ,INTENT(IN), TARGET    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KGMV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,KGMVS) 
REAL(KIND=JPRB)    ,INTENT(IN), TARGET    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KGFL) 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KGT1
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KAUX
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PAUX(YDGEOMETRY%YRDIM%NPROMA,KAUX)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PGPP(YDGEOMETRY%YRDIM%NPROMA,KGT1)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JROF, JLEV, IBL

REAL(KIND=JPRB) :: ZDTDZSG    ! dT/g.dZ

REAL(KIND=JPRB) :: ZTSI(YDGEOMETRY%YRDIM%NPROMA)           ! Interpolated surface temperature
REAL(KIND=JPRB), POINTER :: ZXTCLS(:), ZXRHCLS(:), ZXQCLS(:)
REAL(KIND=JPRB), POINTER :: ZXUCLS(:), ZXVCLS(:), ZXNUCLS(:), ZXNVCLS(:)
REAL(KIND=JPRB), POINTER :: ZGRHL(:,:)
REAL(KIND=JPRB), POINTER :: ZT(:,:)
REAL(KIND=JPRB), POINTER :: ZTL(:,:)
REAL(KIND=JPRB), POINTER :: ZTM(:,:)
REAL(KIND=JPRB), TARGET  :: ZXDUMM(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPRB),ALLOCATABLE :: ZIPH(:,:), ZIPF(:,:), ZTSEA(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "ctstar.intfb.h"
#include "gphpre.intfb.h"
#include "phymfpos.intfb.h"
#include "pos.intfb.h"
#include "gnhqe_conv_tempe.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VPOS',0,ZHOOK_HANDLE)
ASSOCIATE(YDTFP=>YDAFN%TFP,YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDVAB=>YDGEOMETRY%YRVAB, YDPHY=>YDML_PHY_MF%YRPHY, &
 & YDDIMF=>YDML_GCONF%YRDIMF,YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LMPHYS=>YDPHY%LMPHYS, &
 & MXTCLS=>YDXFU%YXFUPT%MXTCLS, MXRHCLS=>YDXFU%YXFUPT%MXRHCLS, MXQCLS=>YDXFU%YXFUPT%MXQCLS,&
 & MXUCLS=>YDXFU%YXFUPT%MXUCLS, MXVCLS=>YDXFU%YXFUPT%MXVCLS, MXNUCLS=>YDXFU%YXFUPT%MXNUCLS, MXNVCLS=>YDXFU%YXFUPT%MXNVCLS, &
 & SD_DI=>YDSURF%SD_DI, SD_VD=>YDSURF%SD_VD, SD_VF=>YDSURF%SD_VF, &
 & SD_VV=>YDSURF%SD_VV, SD_XA=>YDSURF%SD_XA, SD_XR=>YDSURF%SD_XR, &
 & SP_RR=>YDSURF%SP_RR, SP_SG=>YDSURF%SP_SG, YSD_VD=>YDSURF%YSD_VD, &
 & YSD_VF=>YDSURF%YSD_VF, YSP_RR=>YDSURF%YSP_RR)
!     ------------------------------------------------------------------

!*       1.    PREPARATIONS
!              ------------

IBL=(KSTGLO-1)/NPROMA+1

ZDTDZSG=RDTDZ1/RG

IF (YDDYNA%LNHQE) THEN
  ! Conversion Tt -> T and grad(Tt) -> grad (T) with some approximations.
  ! Below calculations, in particuliar POS and PHYMFPOS, are assumed to work with T and grad(T),
  ! and must ignore the 'NHQE' aspects. POS and PHYMFPOS work like with the NHEE model.
  ALLOCATE(ZT(NPROMA,NFLEVG))
  ALLOCATE(ZTL(NPROMA,NFLEVG))
  ALLOCATE(ZTM(NPROMA,NFLEVG))
  DO JLEV=1,NFLEVG
    ZT (:,JLEV) = PGMV(:,JLEV,YDIN_GMV%MT)
    ZTL(:,JLEV) = PGMV(:,JLEV,YDIN_GMV%MTL)
    ZTM(:,JLEV) = PGMV(:,JLEV,YDIN_GMV%MTM)
  ENDDO
  CALL GNHQE_CONV_TEMPE(YDGEOMETRY,YDDYNA%NPDVAR,.TRUE.,1,KST,KEND,PGMV(1,1,YDIN_GMV%MSPD),PGMVS(1,YDIN_GMV%MSP),ZT,&
   & KDDER=2,PQCHAL=PGMV(1,1,YDIN_GMV%MSPDL),PQCHAM=PGMV(1,1,YDIN_GMV%MSPDM),PTL=ZTL,PTM=ZTM)
ELSE
  ZT  => PGMV(:,1:NFLEVG,YDIN_GMV%MT)
  ZTL => PGMV(:,1:NFLEVG,YDIN_GMV%MTL)
  ZTM => PGMV(:,1:NFLEVG,YDIN_GMV%MTM)
ENDIF

! ky: from this line, one should use (ZT,ZTL,ZTM), not the content of PGMV.

IF (LMPHYS) THEN
  DO JROF=KST,KEND
    ZTSI(JROF)=SP_RR(JROF,YSP_RR%YT%MP0,IBL)
  ENDDO
ELSE
  ! ECMWF: only one T profile
  !  ZTSI extrapolated from T(NLEXTRAP)
  ! Compute surface temperatures :
  ALLOCATE(ZIPH(NPROMA,0:NFLEVG))
  ALLOCATE(ZIPF(NPROMA,NFLEVG))
  ALLOCATE(ZTSEA(NPROMA))
  DO JROF=KST,KEND
    ZIPH(JROF,NFLEVG)=EXP(PGMVS(JROF,YDIN_GMV%MSP))
  ENDDO
  CALL GPHPRE(NPROMA,NFLEVG,KST,KEND,YDVAB,YDCVER,ZIPH,PRESF=ZIPF)
  CALL CTSTAR(NPROMA,KST,KEND,ZT(:,NLEXTRAP),ZIPH(1,NFLEVG),&
   & ZIPF(1,NLEXTRAP),YDGEOMETRY%YROROG(IBL)%OROG(:),ZTSI,ZTSEA)  
  DEALLOCATE(ZIPH)
  DEALLOCATE(ZIPF)
  DEALLOCATE(ZTSEA)
ENDIF

!*       2.    PERFORM VERTICAL POST-PROCESSING
!              --------------------------------

CALL POS(YDQTYPE,YDNAMFPSCI,YDTFP,LDHPOS,YDFPVAB,YDGEOMETRY,YDSURF,YDML_GCONF,YDDYN,YDDYNA,YDPHY,NPROMA,KST,KEND,&
 & KPXLEV,KGFL,KGT1,KAUX,NFLEVG,JPNPPM,CDCONF,&
 & YDIN_GFL,YDGFL,YDGEOMETRY%YRGSGEOM(IBL)%RCORI(:),YDGEOMETRY%YRGSGEOM(IBL)%GM(:),YDCSGEOM(IBL)%RATATH(:), &
 & YDGEOMETRY%YROROG(IBL)%OROG(:),YDGEOMETRY%YROROG(IBL)%OROGL(:),YDGEOMETRY%YROROG(IBL)%OROGM(:),&
 & PXLEV,ZTSI,PGMV(1,1,YDIN_GMV%MEDOT),&
 & PGMV(1,1,YDIN_GMV%MU),PGMV(1,1,YDIN_GMV%MUL),PGMV(1,1,YDIN_GMV%MDIV),&
 & PGMV(1,1,YDIN_GMV%MV),PGMV(1,1,YDIN_GMV%MVL),PGMV(1,1,YDIN_GMV%MVOR),&
 & ZT,ZTL,ZTM,&
 & PGMV(1,1,YDIN_GMV%MSPD),PGMV(1,1,YDIN_GMV%MSPDL),PGMV(1,1,YDIN_GMV%MSPDM),&
 & PGMV(1,1,YDIN_GMV%MSVD),PGMV(1,1,YDIN_GMV%MNHX),&
 & PGMVS(1,YDIN_GMV%MSP),PGMVS(1,YDIN_GMV%MSPL),PGMVS(1,YDIN_GMV%MSPM),PGFL,&
 & SD_XA(:,:,:,IBL),SD_XR(:,:,:,IBL),SD_DI(:,:,:,IBL),KCUFNR,PMCUF(:,:),&
 & PGPP,PAUX)

!*       3.    PERFORM PHYSICO-DYNAMIC POST-PROCESSING OVER THE MODEL SURFACE
!              --------------------------------------------------------------

IF(YDGFL%LLG) THEN
  IF (YDGFL%LLH) THEN
    ALLOCATE(ZGRHL(NPROMA,NFLEVG))
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZGRHL(JROF,JLEV)=PGFL(JROF,JLEV,YDIN_GFL%IG)+PGFL(JROF,JLEV,YDIN_GFL%IH)
      ENDDO
    ENDDO
  ELSE
    ZGRHL => PGFL(:,1:NFLEVG,YDIN_GFL%IG)
  ENDIF
ELSE
  ALLOCATE(ZGRHL(NPROMA,NFLEVG))
  ZGRHL(:,:)=0._JPRB
ENDIF

ZXDUMM(:)=HUGE(1._JPRB)
IF (MXTCLS > 0) THEN
  ZXTCLS => YDXFU%XFUBUF(:,MXTCLS,IBL)
ELSE
  ZXTCLS => ZXDUMM(:)
ENDIF
IF (MXRHCLS > 0) THEN
  ZXRHCLS => YDXFU%XFUBUF(:,MXRHCLS,IBL)
ELSE
  ZXRHCLS => ZXDUMM(:)
ENDIF
IF (MXQCLS > 0) THEN
  ZXQCLS => YDXFU%XFUBUF(:,MXQCLS,IBL)
ELSE
  ZXQCLS => ZXDUMM(:)
ENDIF
IF (MXUCLS > 0) THEN
  ZXUCLS => YDXFU%XFUBUF(:,MXUCLS,IBL)
ELSE
  ZXUCLS => ZXDUMM(:)
ENDIF
IF (MXVCLS > 0) THEN
  ZXVCLS => YDXFU%XFUBUF(:,MXVCLS,IBL)
ELSE
  ZXVCLS => ZXDUMM(:)
ENDIF
IF (MXNUCLS > 0) THEN
  ZXNUCLS => YDXFU%XFUBUF(:,MXNUCLS,IBL)
ELSE
  ZXNUCLS => ZXDUMM(:)
ENDIF
IF (MXNVCLS > 0) THEN
  ZXNVCLS => YDXFU%XFUBUF(:,MXNVCLS,IBL)
ELSE
  ZXNVCLS => ZXDUMM(:)
ENDIF

CALL PHYMFPOS(YDQTYPE,YDNAMFPSCI,YDAFN,LDHPOS,YDGEOMETRY,YDDYNA,YDSURF,YDML_GCONF%YRRIP,YDML_PHY_MF,NPROMA,KST,KEND,KAUX, &
 & KGT1,KGFL,YDIN_GFL,YDGFL,YDGEOMETRY%YROROG(IBL)%OROG(:),YDGEOMETRY%YRGSGEOM(IBL)%GM(:),YDGEOMETRY%YRGSGEOM(IBL)%GEMU(:), &
 & YDGEOMETRY%YRGSGEOM(IBL)%GELAM(:),PGMV(1,1,YDIN_GMV%MU),PGMV(1,1,YDIN_GMV%MV),&
 & PGMV(1,1,YDIN_GMV%MDIV),&
 & ZT,PGMV(1,1,YDIN_GMV%MSPD),ZGRHL,PGFL,PGMVS(1,YDIN_GMV%MSP),&
 & SP_RR(:,:,IBL),SP_SG(:,:,:,IBL),SD_VF(:,:,IBL),SD_VV(:,:,IBL),&
 & ZXTCLS,ZXRHCLS,ZXQCLS,ZXUCLS,ZXVCLS,ZXNUCLS,ZXNVCLS,PAUX,PGPP)

IF (YDGFL%LLH .OR. .NOT. YDGFL%LLG) THEN
  ! ZGRHL will have been ALLOCATEd locally
  DEALLOCATE(ZGRHL)
  NULLIFY(ZGRHL)
! ELSE ZGRHL will have been pointed into PGFL so do not deallocate.
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VPOS',1,ZHOOK_HANDLE)

END SUBROUTINE VPOS
