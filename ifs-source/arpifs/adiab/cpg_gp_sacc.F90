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

#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_GP_SACC(YDGEOMETRY,YDMODEL,YDGMV,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & KST,KEND,KSTGLO,LDGW,LDLDIAB,LDMPA,&
 & KIBL,POROG,POROGL,POROGM,PRE0L,PRE0M,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL,PGMV,PRE0,PRE9,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & PRE0F,PXYB0,PUVH0,&
 & PHI0,PHIF0,PHI0FL,PHI0FM,PRCP0,PCTY0,PGWFT0,PGWHT0,PKENE0,&
 & PRT0L,PRT0M,&
 & PRE9F,PXYB9,PHI9,PHIF9,PRCP9,PGWFT9,PGWHT9,&
 & PATND,&
 & PNHXT0,PNHXT9)

!**** *CPG_GP_SACC* - Grid point calculations at instants t and t-dt for shallow atmosphere
!****                 complete-Coriolis equations of Tort & Dubos (2013).

!     Purpose.
!     --------
!      Grid-point calculations at instants t and t-dt for shallow atmosphere complete-coriolis model.
!      Computes some intermediate quantities.
!      Computes some Lagrangian equation RHS.

!**   Interface.
!     ----------
!        *CALL* *CPG_GP_SACC(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTGLO    : global offset.
!        LDGW      : compute [gw].
!        LDLDIAB   : .T. if complete physics is activated and predictor step.
!        LDMPA     : AROME upper-air physics.
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY
!        POROG     : surface orography.
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0L     : zonal component of "grad prehyds" at t.
!        PRE0M     : meridian component of "grad prehyds" at t.

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : unified_treatment grid-point fields at t
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at t-dt.

!     OUTPUT:
!     -------
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PHI0FL    : zonal component of "grad (gz)" at t at full levels.
!        PHI0FM    : meridian component of "grad (gz)" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PGWFT0    : [gw] at full levels at t.
!        PGWHT0    : [gw] at half levels at t.
!        PKENE0    : kinetic energy at t.
!        PRT0L     : zonal component of "grad RT" at full levels.
!        PRT0M     : meridian component of "grad RT" at full levels.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PXYB9     : contains pressure depth, "delta", "alpha" at t-dt.
!        PHI9      : geopotential height "gz" at t-dt at half levels.
!        PHIF9     : geopotential height "gz" at t-dt at full levels.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGWFT9    : [gw] at full levels at t-dt.
!        PGWHT9    : [gw] at half levels at t-dt.
!        PATND     : adiabatic Lagrangian tendencies.
!        PNHXT0    : term 'NHX' at t.
!        PNHXT9    : term 'NHX' at t-dt.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     ----------
!     Author.
!     -------
!      I. Polichtchouk (2021): Adepted CP_GP_HYD to follow LSACC formulation of
!                              the hydrostatic shallow-atm complete-coriolis model of Tort & Dubos (2013)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV   , ONLY : TGMV
USE MODEL_MOD    , ONLY : MODEL

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LSFORC
USE YOMCST   , ONLY : RD, RV
USE YOMCT3   , ONLY : NSTEP
USE YOMLSFORC, ONLY : LSW_FRC, LSOMEGA_FRC
USE YOMGWDIAG, ONLY : UPDATE_GWDIAG,LGWDIAGS_ON
USE INTDYN_MOD,ONLY : YYTHW0, YYTHW9,&
 & YYTCTY0, YYTXYBDER0, YYTRCP0, YYTRCP9, YYTXYB0, YYTXYB9

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
LOGICAL           ,INTENT(IN)    :: LDGW
LOGICAL           ,INTENT(IN)    :: LDLDIAB
LOGICAL           ,INTENT(IN)    :: LDMPA
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE0L(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE0M(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVH0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI0FL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI0FM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWFT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWHT0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKENE0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP9%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWFT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWHT9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YRDYNA%YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV,JROF
LOGICAL :: LLDER,LLGWF0,LLGWF9,LLGDWI

!    ==== MISCELLANEOUS ======================
REAL(KIND=JPRB) :: ZXYBDER0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER0%NDIM)
REAL(KIND=JPRB) :: ZRPRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRDT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZR0T9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRDT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDVER0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDVER9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUS0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZVS0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZUS9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZVS9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZPSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGPHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGPHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUVH9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW9%NDIM)
REAL(KIND=JPRB) :: ZDPHYCTY(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZDUNI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)     ! useless input dummy arg

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gphpre.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gprcp.intfb.h"
#include "gprt.intfb.h"
#include "gpcty.intfb.h"
#include "gpcty_forc.intfb.h"
#include "gpgrgeo.intfb.h"
#include "gpgeo_sacc.intfb.h"
#include "gpgrp_sacc.intfb.h"
#include "gphlwi.intfb.h"
#include "gphluv.intfb.h"
#include "gpuvs.intfb.h"
#include "gpxx.intfb.h"
#include "gpgw.intfb.h"
#include "gp_tndlagadiab_uv_sacc.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_GP_SACC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), &
 & YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE, &
 & YDCVER=>YDGEOMETRY%YRCVER, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,NFLEVG=>YDDIMV%NFLEVG, &
 & NDIM=>YGFL%NDIM, YQ=>YGFL%YQ,YRSPEC=>YGFL%YRSPEC, YFORC=>YGFL%YFORC, &
 & NCURRENT_ITER=>YDMODEL%YRML_DYN%YRDYN%NCURRENT_ITER, &
 & NDIMGMV=>YDGMV%NDIMGMV,YT0=>YDGMV%YT0,YT9=>YDGMV%YT9,&
 & LSIMPH=>YDMODEL%YRML_PHY_MF%YRSIMPHL%LSIMPH)

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

! Fill input useless argument ZDUNI with absurd values.
ZDUNI(KST:KEND,1:NFLEVG)=-999999._JPRB

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES AND TENDENCIES.
!              ----------------------------------------------

!---------------------------------------------------
!      3.1         TIME t0 CALCULATIONS
!---------------------------------------------------

!*     3.1.1 COMPUTE PRE0, PXYB0, PRE0F.

CALL GPHPRE(NPROMA,NFLEVG,KST,KEND,YDVAB,YDCVER,PRE0,PXYB=PXYB0,PRESF=PRE0F)

!*     3.1.2 COMPUTE ZRPRE0F AND ZXYBDER0.

IF(YDCVER%LVERTFE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZRPRE0F(JROF,JLEV)=1.0_JPRB/PRE0F(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF
CALL GPGRXYB(NPROMA,KST,KEND,NFLEVG,.FALSE.,YDVAB,YDCVER,PRE0L,PRE0M,PXYB0,ZXYBDER0)

!*     3.1.3 COMPUTE "R", "Cp" AND "kap=R/Cp".

CALL GPRCP(NPROMA,KST,KEND,NFLEVG,PGFL=PGFL,PCP=PRCP0(1,1,YYTRCP0%M_CP),&
 & PR=PRCP0(1,1,YYTRCP0%M_R),PKAP=PRCP0(1,1,YYTRCP0%M_KAP))

!*     3.1.4 COMPUTES "RT" AND ITS GRADIENT.

IF(YRSPEC%LGP) THEN
  CALL GPRT(YDDYNA%LSPRT,NPROMA,KST,KEND,NFLEVG,RD,RV,&
   & PRCP0(1,1,YYTRCP0%M_R),PGMV(1,1,YT0%MT),PGMV(1,1,YT0%MTL),PGMV(1,1,YT0%MTM),&
   & PGFL(1,1,YQ%MPL),PGFL(1,1,YQ%MPM),&
   & ZRT0,PRT0L,PRT0M,PRL=PGFL(1,1,YRSPEC%MPL),PRM=PGFL(1,1,YRSPEC%MPM))
ELSE
  CALL GPRT(YDDYNA%LSPRT,NPROMA,KST,KEND,NFLEVG,RD,RV,&
   & PRCP0(1,1,YYTRCP0%M_R),PGMV(1,1,YT0%MT),PGMV(1,1,YT0%MTL),PGMV(1,1,YT0%MTM),&
   & PGFL(1,1,YQ%MPL),PGFL(1,1,YQ%MPM),&
   & ZRT0,PRT0L,PRT0M)
ENDIF

! ZRDT0 is the (RT) at "t" used in definition of "dver".
IF (YDDYNA%L_RDRY_VD) THEN
  ZRDT0(KST:KEND,1:NFLEVG)=RD*PGMV(KST:KEND,1:NFLEVG,YT0%MT)
ELSE
  ZRDT0(KST:KEND,1:NFLEVG)=ZRT0(KST:KEND,1:NFLEVG)
ENDIF

!*     3.1.7c COMPUTATION OF SOME VERTICAL VELOCITIES.
!*           Solve continuity equation
!            Computation of the vertical velocities "etapt (d prehyd / d eta)"
!            "omega/prehyd" and "W".

!* Mass tendencies from physics as RHS of the continuity equation.
! Specific mass tendencies (unit: s-1) are computed 
! in the last physics call (previous time step in IFS).
! They are stored in GFL YPHYCTY.
IF (YGFL%YPHYCTY%LACTIVE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZDPHYCTY(JROF,JLEV)=PGFL(JROF,JLEV,YGFL%YPHYCTY%MP)
    ENDDO
  ENDDO
ELSE
  ZDPHYCTY(:,:)=0.0_JPRB
ENDIF

CALL GPCTY(YDVFE,YDCVER,NPROMA,KST,KEND,NFLEVG,YDDYNA%LRUBC,YDVAB,YDVETA,&
 & PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PGMV(1,1,YT0%MDIV),PGMV(1,1,YT0%MEDOT),&
 & PXYB0,PRE0L,PRE0M,ZRPRE0F,PCTY0,PDPHYCTY=ZDPHYCTY)

! Diagnostics of gravity wave noise reflected in surface pressure tendency
IF(LGWDIAGS_ON) THEN
  CALL UPDATE_GWDIAG(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,NPROMA,NFLEVG,KEND,NSTEP,KSTGLO,PRE0,PCTY0(1,1,YYTCTY0%M_DIVDP))
ENDIF

! Overwrite PCTY0 for "EVEL" and "VVEL" with large scale forced values in case
! of 1D model (SCUM)
! (in SCUM, these are 0 in GPCTY)
IF ( LSFORC .AND. (LSW_FRC.OR.LSOMEGA_FRC) ) THEN
  CALL GPCTY_FORC(YDGEOMETRY,YDMODEL%YRML_GCONF,YDMODEL%YRML_PHY_MF%YRPHY2,KST,KEND, &
   & PGFL(1,1,YFORC(1)%MP),PRE0F,PGMV(1,1,YT0%MT),PRCP0(1,1,YYTRCP0%M_R),PCTY0(1,0,YYTCTY0%M_EVEL),PCTY0(1,1,YYTCTY0%M_VVEL))
ENDIF

!*     3.1.8 COMPUTES THE GEOPOTENTIAL HEIGHT "gz" AND ITS GRADIENT.

! * "gz" at full levels and half levels.

PHI0(KST:KEND,NFLEVG)=POROG(KST:KEND)
! * calculate "gz" from modified vertical momentum equation of Tort & Dubos (2013)

CALL GPGEO_SACC(NPROMA,KST,KEND,NFLEVG,PHI0,PHIF0,PGMV(1,1,YT0%MT),PGMV(1,1,YT0%MU),&
 & PRCP0(1,1,YYTRCP0%M_R),PXYB0(1,1,YYTXYB0%M_LNPR),PXYB0(1,1,YYTXYB0%M_ALPH),&
 & YDGEOMETRY%YRVERT_GEOM)  

! * "grad gz" at full levels and half levels.

CALL GPGRGEO(YDGEOMETRY,NPROMA,KST,KEND,NFLEVG,&
 & ZRT0,PRT0L,PRT0M,&
 & PXYB0(1,1,YYTXYB0%M_LNPR),PXYB0(1,1,YYTXYB0%M_ALPH),ZXYBDER0,&
 & POROGL,POROGM,&
 & PHI0FL,PHI0FM,ZGPHL,ZGPHM)

!*     3.1.9 COMPUTES KINETIC ENERGY.

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    PKENE0(JROF,JLEV)=0.5_JPRB*(PGMV(JROF,JLEV,YT0%MU)**2+PGMV(JROF,JLEV,YT0%MV)**2)  
  ENDDO
ENDDO

!*     3.1.11 COMPUTES THE PRESSURE FORCE FOR THE RHS OF MOMENTUM EQN.
!*            AS IN TORT & DUBOS (2013)

CALL GPGRP_SACC(YDGEOMETRY,KST,KEND,&
 & PGMV(1,1,YT0%MU),&
 & ZRT0,PRT0L,PRT0M,PRE0L,PRE0M,PXYB0,ZXYBDER0,&
 & ZGPHL,ZGPHM,PHI0FL,PHI0FM,&
 & ZPSGRTL,ZPSGRTM)

! store for next time-step
IF( YDDYNA%LGRADSP ) THEN

  ! adjust field with filter result
  ZSGRTL(KST:KEND,1:NFLEVG)=ZPSGRTL(KST:KEND,1:NFLEVG)&
   & - (PGMV(KST:KEND,1:NFLEVG,YT9%MSGRTL) - PGMV(KST:KEND,1:NFLEVG,YT0%MSGRTL))
  ZSGRTM(KST:KEND,1:NFLEVG)=ZPSGRTM(KST:KEND,1:NFLEVG)&
   & - (PGMV(KST:KEND,1:NFLEVG,YT9%MSGRTM) - PGMV(KST:KEND,1:NFLEVG,YT0%MSGRTM))

  ! store for next time-step unfiltered fields
  PGMV(KST:KEND,1:NFLEVG,YT0%MSGRTL)=ZPSGRTL(KST:KEND,1:NFLEVG)
  PGMV(KST:KEND,1:NFLEVG,YT0%MSGRTM)=ZPSGRTM(KST:KEND,1:NFLEVG)

  ! set new values for current timestep
  ZPSGRTL(KST:KEND,1:NFLEVG)=ZSGRTL(KST:KEND,1:NFLEVG)
  ZPSGRTM(KST:KEND,1:NFLEVG)=ZSGRTM(KST:KEND,1:NFLEVG)

ENDIF

!*     3.1.12 COMPUTES TERM "NHX".

IF (LDGW) THEN
!*     3.1.13 COMPUTES HALF-LEVEL WINDS.

  CALL GPHLWI(YDGEOMETRY%YRDIMV,NPROMA,KST,KEND,PXYB0(1,1,YYTXYB0%M_LNPR),PXYB0(1,1,YYTXYB0%M_ALPH),PUVH0(1,0,YYTHW0%M_WWI))
  CALL GPHLUV(YDGEOMETRY%YRDIMV,NPROMA,KST,KEND,PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PUVH0)
  LLDER=.FALSE.
  CALL GPUVS(NFLEVG,NPROMA,KST,KEND,LLDER,PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),ZUS0,ZVS0)

!*     3.1.13 COMPUTES TERM "NHX".
  CALL GPXX(YDGEOMETRY,NFLEVG,NPROMA,KST,KEND,ZGPHL,ZGPHM,PHI0FL,PHI0FM,&
   & PXYB0(1,1,YYTXYB0%M_LNPR),ZRT0,&
   & PGMV(1,1,YT0%MU),PGMV(1,1,YT0%MV),PUVH0(1,0,YYTHW0%M_UH),PUVH0(1,0,YYTHW0%M_VH),&
   & PNHXT0)
!*     3.1.14 COMPUTES dver.

  IF (YDDYNA%L_RDRY_VD) THEN
    ! dver_hyd=-(R/Rd)*(D+X+(cv/cp)(omega/prehyd))
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZDVER0(JROF,JLEV)=-(PRCP0(JROF,JLEV,YYTRCP0%M_R)/RD)*(PGMV(JROF,JLEV,YT0%MDIV)&
         & +PNHXT0(JROF,JLEV)+(1.0_JPRB-PRCP0(JROF,JLEV,YYTRCP0%M_KAP))*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL))
      ENDDO
    ENDDO
  ELSE
    ! dver_hyd=-(D+X+(cv/cp)(omega/prehyd))
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZDVER0(JROF,JLEV)=-(PGMV(JROF,JLEV,YT0%MDIV)&
         & +PNHXT0(JROF,JLEV)+(1.0_JPRB-PRCP0(JROF,JLEV,YYTRCP0%M_KAP))*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL))
      ENDDO
    ENDDO
  ENDIF
!*     3.1.16 COMPUTES term [gw].

! Caution: must use consistent definitions of RT in ZRDT0 and ZDVER0.
  LLGWF0=.TRUE.
  LLGDWI=.FALSE.

  CALL GPGW(YDGEOMETRY,YDDYNA%LNHDYN,NFLEVG,NPROMA,KST,KEND,LLGWF0,LLGDWI,&
   & POROGL,POROGM,PXYB0(1,1,YYTXYB0%M_LNPR),PXYB0(1,1,YYTXYB0%M_ALPH),ZUS0,ZVS0,&
   & ZRDT0,ZDVER0,PGWHT0,PGWFT0)
ENDIF

!*     3.1.20a COMPUTES [D V/Dt].

! Pressure gradient term for shallow-atmosphere complete-coriolis eqns + explicit Coriolis + Rayleigh friction
CALL GP_TNDLAGADIAB_UV_SACC(YDGEOMETRY,YDGMV,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_DYN%YRDYN,&
 & KST,KEND,YDGSGEOM%RCORI,YDGSGEOM%GEMU,&
 & YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
 & ZPSGRTL,ZPSGRTM,&
 & PGWFT0,PHIF0,&
 & PGMV,PATND(1,1,YDDYNA%YYTTND%M_TNDU),PATND(1,1,YDDYNA%YYTTND%M_TNDV),&
 & PATND(1,1,YDDYNA%YYTTND%M_TNDU_NOC),PATND(1,1,YDDYNA%YYTTND%M_TNDV_NOC))

!*     3.1.24 COMPUTES [DT/Dt].

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDT)=&
     & PRCP0(JROF,JLEV,YYTRCP0%M_KAP)*PGMV(JROF,JLEV,YT0%MT)*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL)
  ENDDO
ENDDO

!---------------------------------------------------
!      3.2         TIME t9 CALCULATIONS
!---------------------------------------------------

IF (.NOT.YDDYNA%LTWOTL .AND. NCURRENT_ITER == 0) THEN

!*     3.2.1 COMPUTE PRE9, PXYB9, PRE9F.

  CALL GPHPRE(NPROMA,NFLEVG,KST,KEND,YDVAB,YDCVER,PRE9,PXYB=PXYB9,PRESF=PRE9F)

!*     3.2.3 COMPUTE "R", "Cp" AND "Kap=R/Cp".

  CALL GPRCP(NPROMA,KST,KEND,NFLEVG,PGFL=PGFL,KGFLTYP=9,PCP=PRCP9(1,1,YYTRCP9%M_CP),&
   & PR=PRCP9(1,1,YYTRCP9%M_R))

!*     3.2.4 COMPUTE "RT".

  ! ZR0T9 is used in calculation of NHX(t-dt).
  ZR0T9(KST:KEND,1:NFLEVG)=PGMV(KST:KEND,1:NFLEVG,YT9%MT)*PRCP0(KST:KEND,1:NFLEVG,YYTRCP0%M_R)

  ! ZRDT9 is the (RT) at "t-dt" used in the definition of "dver":
  IF (YDDYNA%L_RDRY_VD) THEN
    ZRDT9(KST:KEND,1:NFLEVG)=RD*PGMV(KST:KEND,1:NFLEVG,YT9%MT)
  ELSE
    ! use of R_moist at time t to be consistent with ZR0T9 and ZDVER9.
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZRDT9(JROF,JLEV)=PRCP0(JROF,JLEV,YYTRCP0%M_R)*PGMV(JROF,JLEV,YT9%MT)
      ENDDO
    ENDDO
  ENDIF

!*     3.2.8 COMPUTES THE GEOPOTENTIAL HEIGHT "gz" AND ITS GRADIENT.

  ! * "gz" at full levels and half levels.
  IF(LDLDIAB.OR.LSIMPH) THEN
    PHI9(KST:KEND,NFLEVG)=POROG(KST:KEND)
    CALL GPGEO_SACC(NPROMA,KST,KEND,NFLEVG,PHI9,PHIF9,PGMV(1,1,YT9%MT),PGMV(1,1,YT9%MU),&
     & PRCP9(1,1,YYTRCP9%M_R),PXYB9(1,1,YYTXYB9%M_LNPR),PXYB9(1,1,YYTXYB9%M_ALPH),&
     & YDGEOMETRY%YRVERT_GEOM)  
  ENDIF

!*     3.2.10 COMPUTES HALF-LEVEL WINDS.

  IF (LDGW) THEN
    ! * same remark as for time t half-level winds.
    ZUVH9(KST:KEND,1:NFLEVG,YYTHW9%M_WWI)=PUVH0(KST:KEND,1:NFLEVG,YYTHW0%M_WWI)
    CALL GPHLUV(YDGEOMETRY%YRDIMV,NPROMA,KST,KEND,PGMV(1,1,YT9%MU),PGMV(1,1,YT9%MV),ZUVH9)

    LLDER=.FALSE.
    CALL GPUVS(NFLEVG,NPROMA,KST,KEND,LLDER,PGMV(1,1,YT9%MU),PGMV(1,1,YT9%MV),ZUS9,ZVS9)
  ENDIF

!*     3.2.13 COMPUTES TERM "NHX".

  ! caution: must use "R" at time t.
  IF (LDGW) THEN
    CALL GPXX(YDGEOMETRY,NFLEVG,NPROMA,KST,KEND,ZGPHL,ZGPHM,PHI0FL,PHI0FM,&
     & PXYB9(1,1,YYTXYB9%M_LNPR),ZR0T9,&
     & PGMV(1,1,YT9%MU),PGMV(1,1,YT9%MV),ZUVH9(1,0,YYTHW9%M_UH),ZUVH9(1,0,YYTHW9%M_VH),&
     & PNHXT9)
  ENDIF

!*     3.2.14 COMPUTES dver.

  IF (LDGW) THEN
    IF (YDDYNA%L_RDRY_VD) THEN
      ! dver_hyd=-(R/Rd)*(D+X+(cv/cp)(omega/prehyd))
      !  Caution: must use R, Kap and (omega/prehyd) at time t.
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          ZDVER9(JROF,JLEV)=-(PRCP0(JROF,JLEV,YYTRCP0%M_R)/RD)*(PGMV(JROF,JLEV,YT9%MDIV)&
           & +PNHXT9(JROF,JLEV)+(1.0_JPRB-PRCP0(JROF,JLEV,YYTRCP0%M_KAP))*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL))
        ENDDO
      ENDDO
    ELSE
      ! dver_hyd=-(D+X+(cv/cp)(omega/prehyd))
      !  Caution: must use R, Kap and (omega/prehyd) at time t.
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          ZDVER9(JROF,JLEV)=-(PGMV(JROF,JLEV,YT9%MDIV)&
           & +PNHXT9(JROF,JLEV)+(1.0_JPRB-PRCP0(JROF,JLEV,YYTRCP0%M_KAP))*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL))
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!*     3.2.16 COMPUTES term [gw].
  ! Caution: must use consistent definitions of RT in ZRDT9 and ZDVER9.
  LLGWF9=(.NOT.YDDYNA%LTWOTL).AND.LDMPA
  LLGDWI=.FALSE.
  IF (LDGW) THEN
    CALL GPGW(YDGEOMETRY,YDDYNA%LNHDYN,NFLEVG,NPROMA,KST,KEND,LLGWF9,LLGDWI,&
     & POROGL,POROGM,PXYB9(1,1,YYTXYB9%M_LNPR),PXYB9(1,1,YYTXYB9%M_ALPH),ZUS9,ZVS9,&
     & ZRDT9,ZDVER9,PGWHT9,PGWFT9)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP_SACC',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_GP_SACC
