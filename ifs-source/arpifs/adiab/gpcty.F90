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

!OCL  NOEVAL
SUBROUTINE GPCTY(YDVFE,YDCVER,KPROMA,KD,KF,KFLEV,LDRUBC,YDVAB,YDVETA,&
 & PU,PV,PD,PEVT,PXYB,PSPL,PSPM,PRPREF,PCTY,PDPHYCTY)  

!**** *GPCTY* - Computes vertical velocities.

!     Purpose.
!     --------

!      Computes the vertical velocities
!       "etadot (d prehyd / d eta)" at half levels
!       "(omega / prehyd)" at full levels

!      ----- The following discretisations are valid for lvertfe=.false. -----

!      Omitting the "delta m=1" flux precipitation terms, which will be later
!      added in routine "cpmvvps", the discretised expression of
!      "etadot (d prehyd / d eta)" on the half level number "lbar" is:

!       etadot (d prehyd / d eta) [lbar] =
!        B[lbar] * { sum[k=1 to L]
!        (Delta B)[k] * vec(V)[k] * (grad prehyds) }
!        + B[lbar] * { sum[k=1 to L] (Delta prehyd)[k] * (grad vec(V)[k]) }
!        - sum[k=1 to l] (Delta B)[k] * vec(V)[k] * (grad prehyds)
!        - sum[k=1 to l] (Delta prehyd)[k] * (grad vec(V)[k])
!        + { 1 - B[lbar] } * (etadot (d prehyd / d eta))[top]

!      where:
!       - "vec(V)" is the horizontal wind.
!       - "prehyds" is the surface hydrostatic pressure.
!       - "grad" is the "horizontal gradient" operator:
!         grad X = vnabla X = M vnabla' X

!      Omitting the "delta m=1" flux precipitation terms, which will be later
!      added in routine "cpmvvps", the discretised expression of
!      "(omega / prehyd)" on the full level number "l" is:

!       (omega / prehyd)[l] =
!        vec(V)[l] (grad prehyd / prehyd)[l]
!        - delta[l]/(Delta prehyd)[l] * { sum[k=1 to l-1]
!        (Delta B)[k] * vec(V)[k] * (grad prehyds) }
!        - delta[l]/(Delta prehyd)[l] * { sum[k=1 to l-1]
!        (Delta prehyd)[k] * (grad vec(V)[k]) }
!        - alpha[l]/(Delta prehyd)[l] * (Delta B)[l] * vec(V)[l] * (grad prehyds)
!        - alpha[l] * (grad vec(V)[l])
!        + delta[l]/(Delta prehyd)[l] * (etadot (d prehyd / d eta))[top]

!      This routine stores additional quantities:

!      * vertical integral of divergence without the "lrubc" contribution:
!        for the half level number "lbar" its discretised expression is:

!        psdiv[lbar] = sum[k=1 to l]
!        { (Delta B)[k] * vec(V)[k] * (grad prehyds)
!        + (Delta prehyd)[k] * (grad vec(V)[k]) }

!      * vertical integral of divergence with the "lrubc" contribution:
!        for the half level number "lbar" its discretised expression is:

!        psdvbc[lbar] = sum[k=1 to l]
!        { (Delta B)[k] * vec(V)[k] * (grad prehyds)
!        + (Delta prehyd)[k] * (grad vec(V)[k]) }
!        - (etadot (d prehyd / d eta))[top]

!      * divergence term
!        "pdivdp=grad( vec(V) * (Delta prehyd) )"
!        at full levels: for the full level number "l" its discretised expression is: 

!        pdivdp[l] = 
!        (Delta B)[l] * vec(V)[l] * (grad prehyds)
!        + (Delta prehyd)[l] * (grad vec(V)[l])

!      -----------------------------------------------------------------------

!**   Interface.
!     ----------
!        *CALL* *GPCTY(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          KPROMA       : horizontal dimensioning
!          KD           : start of work
!          KF           : depth of work
!          KFLEV        : number of layers
!          LDRUBC       : upper boundary condition switch
!          YDVAB        : contains information about hybrid vertical coordinate
!          YDVETA       : contains information about hybrid vertical coordinate: "eta"
!          PU           : U component of the wind, at full levels
!          PV           : V component of the wind, at full levels
!          PD           : horizontal divergence of the wind, at full levels
!          PEVT         : "etadot (d prehyd / d eta)" at the top
!          PXYB         : contains pressure depth, "delta", "alpha".
!          PSPL         : zonal component of "grad prehyds"
!          PSPM         : meridian component of "grad prehyds"
!          PRPREF       : inverse of full level pressures.

!        * OUTPUT:
!          PCTY         : contains vertical velocities, vertical integral of divergence.

!        * OPTIONAL INPUT:
!          PDPHYCTY     : mass source/sink from physics (previous dt in LAG physics)

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-02-04

!     Modifications.
!     --------------
!   K. Yessad (Jan 2008): complete (LVERCOR,LVERTFE)=(T,T).
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Dec 2011): use YDVAB.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!     ------------------------------------------------------------------

USE YOMVERT   , ONLY : TVFE, TVAB, TVETA
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCVER   , ONLY : TCVER
USE INTDYN_MOD, ONLY : YYTCTY, YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVFE)        ,INTENT(IN)    :: YDVFE
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KD
INTEGER(KIND=JPIM),INTENT(IN)    :: KF
LOGICAL           ,INTENT(IN)    :: LDRUBC
TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TVETA)       ,INTENT(IN)    :: YDVETA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPREF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTY(KPROMA,0:KFLEV,YYTCTY%NDIM)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PDPHYCTY(KPROMA,KFLEV)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIV(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZPSDIV(KPROMA,KFLEV+1)
REAL(KIND=JPRB) :: ZVP(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZT_EVEL(KPROMA,0:KFLEV),ZT_VVEL(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZT_PSDIV(KPROMA,0:KFLEV),ZT_PSDVBC(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZT_DIVDP(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: JLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPCTY',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!     check for non-supported configurations

IF (LDRUBC.AND.YDCVER%LVERTFE) CALL ABOR1('GPCTY: BAD OPTIONS')

!     ------------------------------------------------------------------

!*    1. Computes "vec(V) * grad prehyds" at full levels.
!     ---------------------------------------------------

DO JLEV=1,KFLEV
  ZVP(KD:KF,JLEV)=PU(KD:KF,JLEV)*PSPL(KD:KF)+PV(KD:KF,JLEV)*PSPM(KD:KF)
ENDDO

!     ------------------------------------------------------------------

!*    2. Sum divergence.
!     ------------------

DO JLEV=1,KFLEV
  ZT_DIVDP(KD:KF,JLEV)=&
   & PD(KD:KF,JLEV)*PXYB(KD:KF,JLEV,YYTXYB%M_DELP)+YDVAB%VDELB(JLEV)*ZVP(KD:KF,JLEV)  
  IF (PRESENT(PDPHYCTY)) THEN
    ZT_DIVDP(KD:KF,JLEV)=ZT_DIVDP(KD:KF,JLEV) &
     & -  PDPHYCTY(KD:KF,JLEV)*PXYB(KD:KF,JLEV,YYTXYB%M_DELP)
  ENDIF
ENDDO

ZT_PSDIV(KD:KF,0)=0.0_JPRB

IF(YDCVER%LVERTFE) THEN
  DO JLEV=1,KFLEV
    ZSDIV(KD:KF,JLEV)=ZT_DIVDP(KD:KF,JLEV)*YDVETA%VFE_RDETAH(JLEV)
  ENDDO
  ZSDIV(KD:KF,0)=0.0_JPRB
  ZSDIV(KD:KF,KFLEV+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',KPROMA,KD,KF,KFLEV,ZSDIV,ZPSDIV)

  DO JLEV=1,KFLEV
    ZT_PSDIV(KD:KF,JLEV)=ZPSDIV(KD:KF,JLEV)
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    ZT_PSDIV(KD:KF,JLEV)=ZT_PSDIV(KD:KF,JLEV-1)+ZT_DIVDP(KD:KF,JLEV)
  ENDDO
ENDIF

IF (LDRUBC) THEN
  IF(YDCVER%LVERTFE) THEN
    ZT_PSDVBC(KD:KF,KFLEV)=ZPSDIV(KD:KF,KFLEV+1)-PEVT(KD:KF)
  ELSE
    DO JLEV=0,KFLEV
      ZT_PSDVBC(KD:KF,JLEV)=ZT_PSDIV(KD:KF,JLEV)-PEVT(KD:KF)
    ENDDO
  ENDIF
ELSE
  IF(YDCVER%LVERTFE) THEN
    ZT_PSDVBC(KD:KF,KFLEV)=ZPSDIV(KD:KF,KFLEV+1)
  ELSE
    ZT_PSDVBC(KD:KF,0:KFLEV)=ZT_PSDIV(KD:KF,0:KFLEV)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*    3. Computes "etadot (d prehyd / d eta)".
!     ----------------------------------------

! "etadot (d prehyd / d eta)" is computed at full levels if LVERTFE=T,
!  at half levels otherwise.

IF(YDCVER%LVERTFE) THEN
  DO JLEV=1,KFLEV
    ZT_EVEL(KD:KF,JLEV)=YDVAB%VBF(JLEV)*ZPSDIV(KD:KF,KFLEV+1)-ZT_PSDIV(KD:KF,JLEV)  
  ENDDO
ELSE
  DO JLEV=1,KFLEV-1
    ZT_EVEL(KD:KF,JLEV)=YDVAB%VBH(JLEV)*ZT_PSDIV(KD:KF,KFLEV)-ZT_PSDIV(KD:KF,JLEV)
  ENDDO
ENDIF

ZT_EVEL(KD:KF,0    )=0.0_JPRB
IF(.NOT.YDCVER%LVERTFE) ZT_EVEL(KD:KF,KFLEV)=0.0_JPRB

IF (LDRUBC) THEN
  IF (YDCVER%LVERTFE) THEN
    CALL ABOR1(' GPCTY: VFE not coded for LRUBC.')
  ELSE
    DO JLEV=0,KFLEV
      ZT_EVEL(KD:KF,JLEV)=ZT_EVEL(KD:KF,JLEV)+(1.0_JPRB-YDVAB%VBH(JLEV))*PEVT(KD:KF)
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*    4. Computes "(omega/prehyd)" at full levels.
!     --------------------------------------------

IF(YDCVER%LVERTFE) THEN

  DO JLEV=1,KFLEV
    ZT_VVEL(KD:KF,JLEV)=(-ZT_PSDIV(KD:KF,JLEV)+YDVAB%VBF(JLEV)*ZVP(KD:KF,JLEV))*PRPREF(KD:KF,JLEV)
  ENDDO

ELSE

  DO JLEV=1,KFLEV
    ZT_VVEL(KD:KF,JLEV)=PXYB(KD:KF,JLEV,YYTXYB%M_RTGR)*ZVP(KD:KF,JLEV)&
     & -PXYB(KD:KF,JLEV,YYTXYB%M_RDELP)*PXYB(KD:KF,JLEV,YYTXYB%M_ALPH)*YDVAB%VDELB(JLEV)*ZVP(KD:KF,JLEV)&
     & -PXYB(KD:KF,JLEV,YYTXYB%M_ALPH)*PD(KD:KF,JLEV)  
  ENDDO
  DO JLEV=1,KFLEV
    ZT_VVEL(KD:KF,JLEV)=ZT_VVEL(KD:KF,JLEV)&
     & -PXYB(KD:KF,JLEV,YYTXYB%M_RDELP)*ZT_PSDVBC(KD:KF,JLEV-1)*PXYB(KD:KF,JLEV,YYTXYB%M_LNPR)  
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*    5. Final memory transfer in PCTY.
!     ---------------------------------

PCTY(KD:KF,0:KFLEV,YYTCTY%M_EVEL)=ZT_EVEL(KD:KF,0:KFLEV)
PCTY(KD:KF,1:KFLEV,YYTCTY%M_VVEL)=ZT_VVEL(KD:KF,1:KFLEV)
PCTY(KD:KF,0:KFLEV,YYTCTY%M_PSDIV)=ZT_PSDIV(KD:KF,0:KFLEV)
PCTY(KD:KF,0:KFLEV,YYTCTY%M_PSDVBC)=ZT_PSDVBC(KD:KF,0:KFLEV)
PCTY(KD:KF,1:KFLEV,YYTCTY%M_DIVDP)=ZT_DIVDP(KD:KF,1:KFLEV)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPCTY',1,ZHOOK_HANDLE)
END SUBROUTINE GPCTY

