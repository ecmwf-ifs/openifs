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
SUBROUTINE GPPVO(YDVAB,KPROMA,KSTART,KPROF,KFLEV,&
 & PRESF,PRDELP,PKAP,PRCORI,&
 & PVOR,PU,PV,PT,PTM,PTL,PSPM,PSPL,PVO,PTETA)

!**** *GPPVO* - COMPUTES POTENTIAL VORTICITY AND POTENTIAL TEMPERATURE.

!     PURPOSE.
!     --------

!      Computes potential temperature and potential vorticity at model full levels.
!      Computes the horizontal derivatives of the potential temperature,
!      and the absolute vorticity as intermediate variables.

!      Potential temperature "pteta" is defined by:
!       pteta = T (prehyd / 100000Pa)**(R/cp)

!      where:
!      - "T" is the temperature
!      - "prehyd" is the hydrostatic pressure
!      - "R" is the air constant (including moisture effects)
!      - "cp" is the constant pressure calorific capacity for moist air.

!      Discretisation at full level "l" of "pteta" follows:
!       pteta[l] = T[l] ( prehyd[l] / 100000Pa )**( R[l] / cp[l] )

!      Discretisation at full level "l" of "grad pteta" follows
!      (horizontal derivatives of (R/cp) are assumed to be negligible):
!       (grad pteta)[l] = pteta[l] * 
!       { (grad T)[l] / T[l]
!       + ( R[l] / cp[l] ) ( grad prehyd / prehyd )[l] }

!      Remark: the current code for ( grad prehyd / prehyd )[l] in part 2
!      gives a discretisation which is not consistent with the one used
!      elsewhere in the code (gradient pressure term, omega/prehyd)
!      which has normally to use the array 'prtgr' computed in 'gpxyb'
!      ( prtgr*(pspl;pspm) yields the good discretisation of vector
!      ( grad prehyd / prehyd ) on layers ). This has to changed in the future.

!      Absolute vorticity "zeta_abs" is defined by:
!       zeta_abs = zeta + f
!      where "zeta" is the (geographic) relative vorticity and "f" is the Coriolis parameter.

!      Potential vorticity "PV" is discretised at full levels, and is defined by:
!       PV = g ( d V / d prehyd ) (grad pteta).vec(i)
!          - g ( d U / d prehyd ) (grad pteta).vec(j)
!          - g zeta_abs ( d pteta / d prehyd )

!      where:
!       - "g" is the gravitation constant.
!       - "U,V" are the components of the horizontal wind. 
!       - "(vec(i),vec(j))" are the unit vectors of the local repere.
!       - "Omega" is the angular velocity of the Earth rotation.
!       - "theta" is the geographic latitude.

!      Vertical derivatives ( d X / d prehyd ) are discretised at full level "l"
!      as follows (X = U; V or "pteta"):
!       ( d X / d prehyd )[l] = 0.5 (X[l+1]-X[l-1])/(Delta prehyd)[l]
!       ( d X / d prehyd )[l=1] = (X[l=2]-X[l=1])/(Delta prehyd)[l=1]
!       ( d X / d prehyd )[l=L] = (X[l=L]-X[l=L-1])/(Delta prehyd)[l=L]
!      This discretisation is always done in finite differences.

!**   INTERFACE.
!     ----------
!        *CALL* *GPPVO(...)

!        EXPLICIT ARGUMENTS :
!        --------------------
!        * INPUT:
!          YDVAB           : TVAB structure (part of geometry).
!          KPROMA          : horizontal dimensioning
!          KSTART to KPROF : depth of work
!          KFLEV           : number of levels
!          PRESF           : hydrostatic pressure "prehyd" at full levels
!          PRDELP          : 1/[Delta prehyd] at full levels
!          PKAP            : "kappa=R/cp" at full levels (previously computed by GPRCP)
!          PRCORI          : Coriolis parameter "f=2 Omega sin(theta)"
!          PVOR            : relative vorticity "zeta" at full levels
!          PU              : U-component of the horizontal wind, at full levels
!          PV              : V-component of the horizontal wind, at full levels
!          PT              : temperature at full levels
!          PTM             : merid grad of temperature "(grad T).vec(j)" at full levels
!          PTL             : zonal grad of temperature "(grad T).vec(i)" at full levels
!          PSPM            : meridian gradient of surface hydrostatic pressure "(grad prehyds).vec(j)"
!          PSPL            : zonal gradient of surface hydrostatic pressure "(grad prehyds).vec(i)"

!        * OUTPUT:
!          PVO             : potential vorticity "PV" at full levels
!          PTETA           : potential temperature "pteta" at full levels

!        IMPLICIT ARGUMENTS :
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!       None.
!       Called by CPG and POS.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!      ORIGINAL : 88-02-04

!     MODIFICATIONS.
!     --------------
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad (Dec 2008): remove useless dummy arguments
!      K. Yessad (March 2017): simplify; remove deep-layer effects
!     ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVAB
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RATM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVO(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTETA(KPROMA,KFLEV) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTETAL(KPROMA,KFLEV),ZTETAM(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPRESL(KPROMA,KFLEV),ZPRESM(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZUZ(KPROMA,KFLEV),ZVZ(KPROMA,KFLEV),ZTETAZ(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZGS2, ZGS2DP, ZGSDPB, ZGSDPT, ZKSPRES, ZUSPT, ZUSRATM, ZVBF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPVO',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    COMPUTES PRESSURE W-E AND N-S DERIVATIVES.
!              ------------------------------------------

DO JLEV=1,KFLEV
  ZVBF=(YDVAB%VBH(JLEV-1)+YDVAB%VBH(JLEV))*0.5_JPRB
  DO JROF=KSTART,KPROF
    ZPRESM(JROF,JLEV) = ZVBF*PSPM(JROF)
    ZPRESL(JROF,JLEV) = ZVBF*PSPL(JROF)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       2. COMPUTES POTENTIAL TEMPERATURE AND ITS W-E AND N-S DERIVATIVE.
!           --------------------------------------------------------------

ZUSRATM=1.0_JPRB/RATM
DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
    ZUSPT  = 1.0_JPRB/PT(JROF,JLEV)
    ZKSPRES= PKAP(JROF,JLEV)/PRESF(JROF,JLEV)
    PTETA(JROF,JLEV)=&
     & +PT(JROF,JLEV)*(PRESF(JROF,JLEV)*ZUSRATM)**(-PKAP(JROF,JLEV))  
    ZTETAM(JROF,JLEV)=PTETA(JROF,JLEV)*&
     & (PTM(JROF,JLEV)*ZUSPT-ZKSPRES*ZPRESM(JROF,JLEV))  
    ZTETAL(JROF,JLEV)=PTETA(JROF,JLEV)*&
     & (PTL(JROF,JLEV)*ZUSPT-ZKSPRES*ZPRESL(JROF,JLEV))  
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTES VERTICAL GRADIENTS.
!              ----------------------------

!     --- LEVELS 2 TO KFLEV-1 ----

ZGS2 = RG*0.5_JPRB

DO JLEV=2,KFLEV-1
  DO JROF=KSTART,KPROF
    ZGS2DP=ZGS2*PRDELP(JROF,JLEV)
    ZUZ   (JROF,JLEV)=(PU   (JROF,JLEV+1)-PU   (JROF,JLEV-1))*ZGS2DP
    ZVZ   (JROF,JLEV)=(PV   (JROF,JLEV+1)-PV   (JROF,JLEV-1))*ZGS2DP
    ZTETAZ(JROF,JLEV)=(PTETA(JROF,JLEV+1)-PTETA(JROF,JLEV-1))*ZGS2DP
  ENDDO
ENDDO

!     --- LEVELS 1 AND KFLEV ----

DO JROF=KSTART,KPROF
  ZGSDPT=RG*PRDELP(JROF,1)
  ZUZ   (JROF,1)=(PU   (JROF,2)-PU   (JROF,1))*ZGSDPT
  ZVZ   (JROF,1)=(PV   (JROF,2)-PV   (JROF,1))*ZGSDPT
  ZTETAZ(JROF,1)=(PTETA(JROF,2)-PTETA(JROF,1))*ZGSDPT
  ZGSDPB=RG*PRDELP(JROF,KFLEV)
  ZUZ   (JROF,KFLEV)=(PU   (JROF,KFLEV)-PU   (JROF,KFLEV-1))*ZGSDPB
  ZVZ   (JROF,KFLEV)=(PV   (JROF,KFLEV)-PV   (JROF,KFLEV-1))*ZGSDPB
  ZTETAZ(JROF,KFLEV)=(PTETA(JROF,KFLEV)-PTETA(JROF,KFLEV-1))*ZGSDPB
ENDDO

!     ------------------------------------------------------------------

!*       4.    COMPUTES POTENTIAL VORTICITY.
!              -----------------------------

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
    PVO(JROF,JLEV)=&
     & -(PVOR(JROF,JLEV)+PRCORI(JROF))*ZTETAZ(JROF,JLEV)&
     & +ZVZ(JROF,JLEV)*ZTETAL(JROF,JLEV)&
     & -ZUZ(JROF,JLEV)*ZTETAM(JROF,JLEV)  
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPVO',1,ZHOOK_HANDLE)
END SUBROUTINE GPPVO
