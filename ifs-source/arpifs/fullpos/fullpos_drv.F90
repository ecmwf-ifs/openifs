! (C) Copyright 1989- Meteo-France.

SUBROUTINE FULLPOS_DRV(YDFPOS,YDGEOMETRY,YDGMV,YDGFL,YDSURF,YDCFU,YDXFU,PCFUBUF,PXFUBUF,KCUFNR,PMCUFGP, &
 & YDMODEL,PTSTEP,KSTEP,KSTOP,KFPLAG,LDOCEDELAY,LDOCE,YDFPDATA,KFPPHY,KMFPPHY)

!**** *FULLPOS_DRV*  - Fullpos process

!     Purpose.
!     --------
!        To process Fullpos

!**   Interface.
!     ----------
!        *CALL* *FULLPOS_DRV(...)

!        Explicit arguments :
!        --------------------
!           YDGEOMETRY : input model geometry
!           PTSTEP : model time step
!           KSTEP : current model time step
!           KSTOP : last model time step
!           KFPLAG : lagged/unlagged mode to call the post-processing (for physics)
!                    = 0 : unlagged mode
!                    = 1 : lagged mode, part 1 (physics + dynamics)
!                    = 2 : lagged mode, part 1 (lagged physics)
!           KFPPHY   : number of requested surface fields (ECMWF only)
!           KMFPPHY  : GRIB codes of requested surface fields (ECMWF only)
!           LDOCEDELAY : Post-processing physics fields at the end of the step
!           LDOCE : Post-processing ocean fields if called in ocean mode


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

!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE FA_MOD       , ONLY : JD_SIZ
USE YOMCT0       , ONLY : LARPEGEF, N3DINI
USE YOMOPH0      , ONLY : LINC
USE YOMFPC       , ONLY : L_READ_MODEL_DATE, LFPMOIS
USE FULLPOS      , ONLY : TFPOS, TFPDATA
USE YOMGMV       , ONLY : TGMV
USE YOMGFL       , ONLY : TGFL
USE YOMCFU       , ONLY : TCFU
USE YOMXFU       , ONLY : JPFUXT, TXFU
USE YOMLUN       , ONLY : NULOUT
USE YOMRIP0      , ONLY : NINDAT
USE YOMARG       , ONLY : NUDATE
USE TYPE_FPFIELDS, ONLY : TFPFIELDS, DEALLOC_FPFIELDS
USE YOMFP4L , ONLY : TRQFP
USE YOM4FPOS, ONLY : TRQFPDYN
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE SPECTRAL_FIELDS_MOD

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(TFPOS)       , TARGET, INTENT(IN) :: YDFPOS
TYPE(GEOMETRY)    ,         INTENT(IN) :: YDGEOMETRY
TYPE(TGMV)        ,         INTENT(IN) :: YDGMV
TYPE(TGFL)        ,         INTENT(IN) :: YDGFL
TYPE(TSURF)       ,         INTENT(IN) :: YDSURF
TYPE(TCFU)        ,         INTENT(IN) :: YDCFU
TYPE(TXFU)        ,         INTENT(IN) :: YDXFU
REAL(KIND=JPRB)   ,         INTENT(IN) :: PCFUBUF(YDGEOMETRY%YRDIM%NPROMA,YDCFU%NFDCFU,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,         INTENT(IN) :: PXFUBUF(YDGEOMETRY%YRDIM%NPROMA,YDXFU%NFDXFU,YDGEOMETRY%YRDIM%NGPBLKS)
INTEGER(KIND=JPIM),         INTENT(IN) :: KCUFNR
REAL(KIND=JPRB)   ,         INTENT(IN) :: PMCUFGP(YDGEOMETRY%YRDIM%NPROMA,KCUFNR,YDGEOMETRY%YRDIM%NGPBLKS)
TYPE(MODEL)       ,         INTENT(IN) :: YDMODEL
REAL(KIND=JPRB)   ,         INTENT(IN) :: PTSTEP
INTEGER(KIND=JPIM),         INTENT(IN) :: KSTEP
INTEGER(KIND=JPIM),         INTENT(IN) :: KSTOP
INTEGER(KIND=JPIM),         INTENT(IN) :: KFPLAG
LOGICAL           ,         INTENT(IN) :: LDOCEDELAY
LOGICAL           ,         INTENT(IN) :: LDOCE
TYPE(TFPDATA)     ,         INTENT(INOUT) :: YDFPDATA
INTEGER(KIND=JPIM),         INTENT(OUT) :: KFPPHY
INTEGER(KIND=JPIM),         INTENT(OUT) :: KMFPPHY(:)

!     ------------------------------------------------------------------
LOGICAL :: LLRQCFU, LLRQXFU, LLAGED, LLPHYS, LLUPDSUW, LLOCE
INTEGER(KIND=JPIM) :: IMMCLI, IDTIME, IMODELTIME, IFPDOM, ITER
INTEGER(KIND=JPIM) :: IDATEF(JD_SIZ)
TYPE (TFPFIELDS)   :: YLFPFIELDS
TYPE(TRQFP) :: YLRQPHY
TYPE(TRQFPDYN) :: YLRQDYN

CHARACTER (LEN = 35) ::  CLINE=' '

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "su4fpos.intfb.h"
#include "sufpdata.intfb.h"
#include "allfpos.intfb.h"
#include "suppdate.intfb.h"

#include "fcttim.func.h"

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FULLPOS_DRV',0,ZHOOK_HANDLE)
ASSOCIATE(YDAFN=>YDFPOS%YAFN, YDFPGEOMETRY=>YDFPOS%YFPGEOMETRY, YDFPOFN=>YDFPOS%YFPIOH%YFPOFN(:),YDFPVAB=>YDFPOS%YFPVAB, &
 & YDFPSTRUCT=>YDFPOS%YFPSTRUCT, YDFPWSTD=>YDFPOS%YFPWSTD, YDFPCNT=>YDFPOS%YFPCNT, YDNAMFPIOS=>YDFPOS%YFPIOH%YNAMFPIOS, &
 & YDNAMFPSCI=>YDFPOS%YNAMFPSCI, YDNAMFPINT=>YDFPOS%YNAMFPINT,YDML_GCONF=>YDMODEL%YRML_GCONF)
ASSOCIATE(YDFPGEO_DEP=>YDFPGEOMETRY%YFPGEO_DEP, YDFPGEO=>YDFPGEOMETRY%YFPGEO, YDFPUSERGEO=>YDFPGEOMETRY%YFPUSERGEO, &
 & LFPOSHOR=>YDFPGEOMETRY%LFPOSHOR, NFPXFLD=>YDNAMFPIOS%NFPXFLD,NFPSURFEX=>YDNAMFPSCI%NFPSURFEX, NFPCLI=>YDNAMFPSCI%NFPCLI)
ASSOCIATE(CFPDOM=>YDFPUSERGEO(:)%CFPDOM)
!      -----------------------------------------------------------


!*       1.    VARIOUS INITIALIZATION
!              ---------------------

IFPDOM=SIZE(YDFPUSERGEO)
LLAGED = (KFPLAG >= 2)

LLPHYS=YDMODEL%YRML_PHY_EC%YREPHY%LEPHYS.OR.YDMODEL%YRML_PHY_MF%YRPHY%LMPHYS

! Time-dependent request
CALL GSTATS(29,0)
WRITE(NULOUT,*) ' == Full-Pos : allocate time-dependent arrays == ',CLINE
IMODELTIME=INT(REAL(KSTEP,KIND=JPRB)*PTSTEP+0.5_JPRB)
CALL SU4FPOS(YDNAMFPSCI,YDAFN,YDFPCNT,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NCSNEC,LLPHYS,IMODELTIME,IFPDOM,CFPDOM,YDFPVAB,YDML_GCONF%YGFL, &
 & YDMODEL%YRML_PHY_RAD%YREAERATM,LLAGED,LDOCEDELAY,LDOCE,YLFPFIELDS,YLRQPHY,YLRQDYN)
IF (.NOT.LARPEGEF) THEN
  KFPPHY=YLFPFIELDS%NFPPHY
  KMFPPHY(1:KFPPHY)=YLFPFIELDS%MFPPHY(1:KFPPHY)
ENDIF
CALL GSTATS(29,1)


!*       2.    INITIALIZE DATE FOR I/Os
!              ------------------------

CALL GSTATS(1776,0)

IF (LARPEGEF .AND. KFPLAG < 2) THEN
  ! Compute output date for FA. For grib the model date is used.
  IF (L_READ_MODEL_DATE) THEN
    IDTIME=0
  ELSE
    IDTIME=1-N3DINI
  ENDIF
  IF (YDXFU%LXFU) THEN
    CALL SUPPDATE(YDML_GCONF%YRRIP,KSTEP,IDTIME,KEVENTTS=YDXFU%NRAZTS,KDIME=JPFUXT, &
 & KFREVENT=YDXFU%NFRRAZ,K1EVENT=YDXFU%N1RAZ,KPPDATE=IDATEF)
  ELSE
    CALL SUPPDATE(YDML_GCONF%YRRIP,KSTEP,KDTIME=IDTIME,KDIME=1,KPPDATE=IDATEF)
  ENDIF
ENDIF

IF (.NOT.LLAGED) THEN

!*       3.    INITIALIZE SURFACE-RELATED AUXILARY DATA
!              ----------------------------------------

  IF (YDFPDATA%LFPUPDCLI) THEN
    IF (LFPMOIS) THEN
    ! Date of the starting file (no need to change the climatology file)
      IMMCLI=NMM(NUDATE)
    ELSE
    ! Date of the model (more accurate, but needs to update the climatology file)
      IMMCLI=NMM(NINDAT)
    ENDIF
  ELSE
    IMMCLI=0
  ENDIF
  ITER=0
  LLUPDSUW=YDFPDATA%LFPUPDSUW
  CALL SUFPDATA(YDFPOS,YDGEOMETRY,YDSURF,YDMODEL,IMMCLI,ITER,YDFPDATA,LLUPDSUW)

ELSE
  ! MC consider passing YDFPDATA to model_step instead
  YDFPDATA%YFPOS=>YDFPOS

ENDIF

!*       5.    POST-PROCESSING
!              ---------------

LLRQCFU=(YDCFU%NFDCFU > 0).AND.(YDCFU%LREACFU.OR.KSTEP > 0)
LLRQXFU=(YDXFU%NFDXFU > 0).AND.(YDXFU%LREAXFU.OR.KSTOP > 0)

LLOCE=LDOCEDELAY.AND.LDOCE
CALL ALLFPOS(YLRQDYN,YDFPDATA,YLRQPHY,LLAGED,LLOCE,LLRQCFU,LLRQXFU,YDGEOMETRY,YDGMV,YDGFL, &
 & YDSURF,YDCFU,YDXFU,YDMODEL,PCFUBUF,PXFUBUF,KCUFNR,PMCUFGP,IDATEF,KSTEP, &
 & LDINC=LINC,KPPSTEP=KSTEP,PTSTEP=PTSTEP)


CALL GSTATS(1776,1)

WRITE(NULOUT,*) ' == Full-Pos : release time-dependent arrays == ',CLINE
CALL DEALLOC_FPFIELDS(YLFPFIELDS)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FULLPOS_DRV',1,ZHOOK_HANDLE)
END SUBROUTINE FULLPOS_DRV
