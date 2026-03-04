! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPACHMT(YDML_PHY_MF,LDZ0H,KPROMA,KSTART,KPROF,KFLEV,&
 & PTT0,PQT0,PUT0,PVT0,&
 & PRESH,PRESF,PCP0,PR0,PGEOPH,PGEOPF,&
 & PLSM,PTS,PGZ0F,PGZ0H,PVEG0,PARG,PD2,PIVEG,&
 & PSAB,PSNS,PWS,PWSI,PHV,PQCLS,PRHCLS,PTCLS,PUCLS,PVCLS,&
 & PNUCLS,PNVCLS,PFCLS,PUGST,PVGST,PFGST,PRPCLS,PQS) 

!     PURPOSE.
!     --------
!        Full-POS interface to ACHMT, to compute pressure, humidity,
!        temperature and wind

!**   INTERFACE.
!     ----------
!       *CALL* *FPACHMT*

!        EXPLICIT ARGUMENTS
!        --------------------
!               INPUT :
!               ------
!        KPROMA    : horizontal dimension
!        KSTART    : start of work
!        KPROF     : depth of work
!        KFLEV     : numer of full levels
!        PTT0      : temperature
!        PQT0      : specific moisture
!        PUT0      : geographical wind - U momentum
!        PVT0      : geographical wind - V momentum
!        PRESH     : half level hydrostatic pressure
!        PRESF     : full level hydrostatic pressure
!        PCP0      : Cp
!        PR0       : R (moist air constant)
!        PGEOPH    : geopotential height at half levels
!        PGEOPF    : geopotential height at full levels
!        PLSM      : land/sea mask
!        PTS       : surface temperature
!        PGZ0F     : gravity * roughness length
!        PGZ0H     : gravity * thermal roughness length
!        PVEG0     : fractional cover by vegetation
!        PARG      : silt percentage
!        PD2       : soil depth
!        PIVEG     : type of vegetation
!        PSAB      : sand percentage 
!        PSNS      : snow cover per unit surface
!        PWS       : surface water content
!        PWSI      : surface layer ice content
!        PHV       : resistance to evapotranspiration

!               OUTPUT:
!               ------
!        PQCLS     : specific moisture at 2m
!        PRHCLS    : relative moisture at 2m
!        PTCLS     : temperature at 2m
!        PUCLS     : U-component of wind at 10m
!        PVCLS     : V-component of wind at 10m
!        PNUCLS    : U-component of neutral wind at 10m
!        PNVCLS    : V-component of neutral wind at 10m
!        PFCLS     : Wind velocity of wind at 10m
!        PUGST     : U-component of gusts
!        PVGST     : V-component of gusts
!        PFGST     : gusts
!        PRPCLS    : pressure
!        PQS       : surface pressure

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE* & WITOLD OWCARZ
!      ORIGINAL   : 2000-03-08 by Witold Owcarz

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 2001-01-05 Gusts
!      J.M. Piriou / J.F Geleyn  : 2001-08-23 New definition and use of the
!                      thermal mixing length.
!      M. Bellus    : 2003-04-24 Change of arguments in ACHMT calling
!                                + removing LSOLV and some cleaning
!      M. Tudor  : 2003-10-21 Change args in ACHMT, reintroduce LSOLV
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      A. Alias  : 2007-10-11 New Output in ACHMT (ZRTI=1/R*T)
!      C. Payan  : 2008/11 Neutral Wind Diagnostic (ACHMT)
!      K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!      K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!      K. Yessad (Aug 2009): move GP.. calculations in caller.
!      R. El Khatib 13-Dec-2012 Fix intent attribute
!      09-2018 R. Brozkova, A. Bucanek: MOCON diagnostics in offline fullpos
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
LOGICAL           ,INTENT(IN)    :: LDZ0H
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP0(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR0(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOPH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOPF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0F(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0H(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEG0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PARG(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD2(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIVEG(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSAB(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHV(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRHCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNUCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNVCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUGST(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVGST(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFGST(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRPCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQS(KPROMA)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZGZ0RLF(KPROMA),ZGZ0(KPROMA),ZGZ0H(KPROMA)
REAL(KIND=JPRB) :: ZCD(KPROMA),ZCDN(KPROMA),ZCH(KPROMA),ZCS(KPROMA)
REAL(KIND=JPRB) :: ZNEIJJ(KPROMA),ZVEG(KPROMA)
REAL(KIND=JPRB) :: ZRS(KPROMA),ZQSATS(KPROMA)
REAL(KIND=JPRB) :: ZNBVNO(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZMRIPP(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZFPLSH(KPROMA,0:KFLEV),ZFPLCH(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZXDROV(KPROMA),ZXHROV(KPROMA)
REAL(KIND=JPRB) :: ZDPHIV(KPROMA),ZDPHIT(KPROMA)
REAL(KIND=JPRB) :: ZNEIJG(KPROMA),ZNEIJV(KPROMA),ZWFC(KPROMA)
REAL(KIND=JPRB) :: ZCDROV(KPROMA),ZCHROV(KPROMA),ZDQSTS(KPROMA),ZGWDCS(KPROMA)
REAL(KIND=JPRB) :: ZHQ(KPROMA),ZHU(KPROMA),ZSTAB(KPROMA)
REAL(KIND=JPRB) :: ZRTI(KPROMA)
REAL(KIND=JPRB) :: ZWPMX(KPROMA),ZWSMX(KPROMA),ZWWILT(KPROMA),ZWSAT(KPROMA)

LOGICAL :: LLCLS, LLHMT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "achmt.intfb.h"
#include "actkehmt.intfb.h"
#include "acsolw.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPACHMT',0,ZHOOK_HANDLE)
ASSOCIATE(HVCLS=>YDML_PHY_MF%YRPHY2%HVCLS, HTCLS=>YDML_PHY_MF%YRPHY2%HTCLS, LRAFTUR=>YDML_PHY_MF%YRPHY2%LRAFTUR, &
 & LSNV=>YDML_PHY_MF%YRPHY%LSNV, LCOEFKSURF=>YDML_PHY_MF%YRPHY%LCOEFKSURF, LSOLV=>YDML_PHY_MF%YRPHY%LSOLV)
!-----------------------------------------------------------------------

!*       0.    RELIABITY OF INPUT ARGUMENTS
!              ----------------------------

IF (LSNV) THEN
  CALL ABOR1('FPACHMT : NOT READY FOR LSNV - ZNEIJG,ZNEIJV')
ENDIF

!-----------------------------------------------------------------------

!        1.    COMPUTATION OF T, Q, U AND V
!              ----------------------------

LLHMT=.TRUE.

IF (LSOLV) THEN
  CALL ACSOLW(YDML_PHY_MF%YRPHY1,KSTART,KPROF,KPROMA,PARG,PD2,PLSM,PIVEG,PSAB,LLHMT,ZWFC,ZWPMX,&
   & ZWSAT,ZWSMX,ZWWILT)
ENDIF

ZDPHIV(KSTART:KPROF)=RG*HVCLS
ZDPHIT(KSTART:KPROF)=RG*HTCLS
LLCLS=.FALSE.
 
IF (LCOEFKSURF) THEN
  CALL ACTKEHMT(KSTART,KPROF,KPROMA,KFLEV,LDZ0H.AND.LSOLV,&
   & PGEOPH,PGEOPF,PRESH,PRESF,PCP0,PQT0,PR0,PTT0,PUT0,PVT0,&
   & ZFPLSH,ZFPLCH,ZDPHIT,ZDPHIV,PGZ0F, PGZ0H, ZGZ0RLF, PHV, PLSM,&
   & ZNEIJG,ZNEIJV, PSNS, PTS, PVEG0, ZWFC, PWS, PWSI,&
   & LLCLS, LLHMT,&
   & ZNBVNO,ZMRIPP,&
   & ZCD, ZCDN, ZCDROV, ZCH, ZCHROV, ZCS, ZDQSTS, ZGWDCS,&
   & ZGZ0, ZGZ0H, ZHQ, ZHU, ZNEIJJ, PQCLS, PQS, ZQSATS, PRHCLS,&
   & ZRS, ZRTI, ZSTAB, PTCLS, PUCLS, PVCLS, PNUCLS, PNVCLS, PRPCLS, ZVEG,&
   & ZXDROV, ZXHROV, PUGST, PVGST)  
ELSE
  CALL ACHMT(KSTART,KPROF,KPROMA,KFLEV,LDZ0H.AND.LSOLV,&
   & PGEOPH,PGEOPF,PRESH,PRESF,PCP0,PQT0,PR0,PTT0,PUT0,PVT0,&
   & ZFPLSH,ZFPLCH,ZDPHIT,ZDPHIV,PGZ0F, PGZ0H, ZGZ0RLF, PHV, PLSM,&
   & ZNEIJG,ZNEIJV, PSNS, PTS, PVEG0, ZWFC, PWS, PWSI,&
   & LLCLS, LLHMT,&
   & ZNBVNO,ZMRIPP,&
   & ZCD, ZCDN, ZCDROV, ZCH, ZCHROV, ZCS, ZDQSTS, ZGWDCS,&
   & ZGZ0, ZGZ0H, ZHQ, ZHU, ZNEIJJ, PQCLS, PQS, ZQSATS, PRHCLS,&
   & ZRS, ZRTI, ZSTAB, PTCLS, PUCLS, PVCLS, PNUCLS, PNVCLS, PRPCLS, ZVEG,&
   & ZXDROV, ZXHROV, PUGST, PVGST)  
ENDIF

IF (.NOT.LRAFTUR) THEN 
  CALL ABOR1('FPACHMT:CALL TO ACCLPH BROKEN\MH')
  ! CALL ACCLPH (KSTART,KPROF,KPROMA,1,KFLEV,PGEOPH,PGEOPF,PR0,PTT0,PUT0,PVT0,&
  !  &ZRS,PTS,ICLPH,ZCLPH,PUGST,PVGST,'ACCLPH ')
  ! else they are computed inside ACHMT
ENDIF

PFCLS(KSTART:KPROF)=SQRT((PUCLS(KSTART:KPROF)**2)+(PVCLS(KSTART:KPROF)**2))
PRHCLS(KSTART:KPROF)=MAX(0.0_JPRB,MIN(1.0_JPRB,PRHCLS(KSTART:KPROF)))
PFGST(KSTART:KPROF)=SQRT((PUGST(KSTART:KPROF)**2)+(PVGST(KSTART:KPROF)**2))

!-------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPACHMT',1,ZHOOK_HANDLE)
END SUBROUTINE FPACHMT
