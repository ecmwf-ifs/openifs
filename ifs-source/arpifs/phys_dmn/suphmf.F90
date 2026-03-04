! (C) Copyright 1989- Meteo-France.

!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHMF(YDGEOMETRY,YDMODEL,KULOUT)

!**** *SUPHMF*   - Calls initialization of commons controlling physics
!                  in the Meteo-France version.

!     Purpose.
!     --------
!           Organise the setup of physical constants for Meteo-France
!           physics package.

!**   Interface.
!     ----------
!        *CALL* *SUPHMF(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        None.

!     Method.
!     -------
!        Irrelevant.

!     Externals.
!     ----------

!     SUPHY0
!     SUPHY1
!     SUPHY2
!     SUPHY3
!     SUTOPH

!     Reference.
!     ----------

!     Author.
!     -------
!      J.-F. Geleyn.
!      Original : 91-06-15

!     Modifications.
!     --------------
!      Modified 01-04-02 R. El Khatib setup for CAPE diagnostic
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 04-11-16 Y. Seity : call suphmnh for AROME physics
!      R. Zaaboul 28-Feb-2006: call suparar, suphmpa and suphmse (ex suphmnh)
!      Y. Seity 06-07-10: nfpsurfex in call suphmse (prepsurfex)  
!      E . Bazile 06-09-01 : call to suphmpa if LCVPPKF for ARPEGE/ALADIN.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. Brozkova 09-2018: Modified call to SUCAPE.
!     ------------------------------------------------------------------

!*       1.    Call routines for specific physics' commons setup.
!              --------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "sucape.intfb.h"
#include "suphy0.intfb.h"
#include "suphy1.intfb.h"
#include "suphy2.intfb.h"
#include "suphy3.intfb.h"
#include "sutoph.intfb.h"
#include "sunorgwd.intfb.h"
#include "val923.intfb.h"
#include "suparar.intfb.h"
#include "suphmpa.intfb.h"
#include "suphmse.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUPHMF',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDLAP=>YDGEOMETRY%YRLAP,  &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
  & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,  &
  & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY)

ASSOCIATE(LMSE=>YDARPHY%LMSE, LMPA=>YDARPHY%LMPA, &
 & LSOLV=>YDPHY%LSOLV, LEDKF=>YDPHY%LEDKF, LCVPPKF=>YDPHY%LCVPPKF)
CALL SUPHY0(YDGEOMETRY%YRGEM,YDMODEL%YRML_PHY_MF,KULOUT)
CALL SUPHY1(YDMODEL%YRML_PHY_MF%YRPHY1,YDMODEL%YRML_PHY_MF%YRVDOZ,KULOUT)
CALL SUPHY2(YDVAB,YDGEOMETRY%YRDIMV,YDMODEL%YRML_PHY_MF%YRPHY2,KULOUT)
CALL SUPHY3(YDPHY,YDMODEL%YRML_PHY_MF%YRPHY3,KULOUT)
CALL SUTOPH(YDSTA,YDGEOMETRY%YRDIMV,YDMODEL%YRML_PHY_EC%YREPHY,YDPHY,YDMODEL%YRML_PHY_MF%YRTOPH,KULOUT)

! SETTING CONSTANTS FOR NON-OROGRAPHIC GWD SCHEME
CALL SUNORGWD(YDSTA,YDGEOMETRY%YRDIMV,YDMODEL%YRML_PHY_MF%YRNORGWD,KULOUT)

CALL VAL923(LSOLV)

CALL SUCAPE(YDPHY,YDVAB,YDDIMV,KULOUT)

! setup for AROME physics and SURFEX
CALL SUPARAR(YDGEOMETRY,YDMODEL%YRML_GCONF%YGFL,YDMODEL%YRML_PHY_MF,KULOUT)
IF (LMPA.OR.LCVPPKF.OR.LEDKF) CALL SUPHMPA(YDGEOMETRY,YDMODEL%YRML_DIAG%YRLDDH,YDMODEL%YRML_GCONF, &
 & YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_PHY_MF,KULOUT)
IF (LMSE) CALL SUPHMSE(YDGEOMETRY,YDMODEL,KULOUT)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHMF',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHMF

