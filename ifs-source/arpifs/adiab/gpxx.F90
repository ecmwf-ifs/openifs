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

SUBROUTINE GPXX(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY, KFLEV,KPROMA,KSTART,KEND,PHIHL,PHIHM,PHIFL,PHIFM,&
 & PLNPR,PRT,&
 & PUF,PVF,PUH,PVH,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PX,&
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & PNHPPI&
 & )

! GPXX - Diagnose NHX-term

! Purpose
! -------
!   Diagnose NHX-term 
!    NHX = (pre/(RT)) grad[gz] (d V / d prehyd)
!   It is better to rewrite NHX as follows for the discretisation:
!    NHX = (pre/prehyd) (1/(RT)) grad[gz] ( d V / d (log(prehyd)) )

!   NHX is discretised as follows, at full levels, for FD discretisation:
!   [X]_[l] = [pre/prehyd]_[l] * [ 1/(R_[l] T_[l] delta_[l]) ] *
!             [ grad[gz]_[lbar-1] (V_[l] - V[lbar-1])
!             + grad[gz]_[lbar] (V_[lbar] - V_[l]) ]
!   VFE discretisation is different.

!   This routine can be used in a NHEE model: PNHPPI should be present in this case.

!   This routine can also be used in a hydrostatic model or in a NHQE model:
!   in this case the ratio (pre/prehyd) is equal to 1 and PNHPPI should not be used.

! Interface
! ---------
!   * INPUT:
!   YDGEOMETRY   : structure containing all geometry
!   KFLEV        : number of levels.
!   KPROMA       : length of work
!   KSTART       : start of work
!   KEND         : end of work
!   PHIHL        : zonal component of "grad[gz]" at half levels.
!   PHIHM        : meridian component of "grad[gz]" at half levels.
!   PHIFL        : zonal component of "grad[gz]" at full levels.
!   PHIFM        : meridian component of "grad[gz]" at full levels.
!   PLNPR        : "delta" at full levels.
!   PRT          : (R*Temperature) at full levels.
!   PUF          : U-wind at full levels.
!   PVF          : V-wind at full levels.
!   PUH          : U-wind at half levels.
!   PVH          : V-wind at half levels.

!   * OUTPUT:
!   PX           : NHX-term at full levels.

!   * OPTIONAL INPUT:
!   PNHPPI       : [pre/prehyd] at full levels (NHEE model).

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   06 Dec 2004 K. Yessad (after GNHX).

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   J. Vivoda (Oct 2013): new options for VFE-NH
!   K. Yessad (June 2017): Introduce NHQE model.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIHL(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIHM(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIFL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIFM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PX(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PNHPPI(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZDUF(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZDVF(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZF(KPROMA,0:KFLEV+1)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "verdisint.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPXX',0,ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

IF( YDGEOMETRY%YRCVER%LVERTFE .AND. YDGEOMETRY%YRCVER%LVFE_X_TERM )THEN

  ! vertical derivatives of wind
  IF (YDGEOMETRY%YRCVER%NVFE_DERBC==0.OR.YDGEOMETRY%YRCVER%NVFE_DERBC==1) THEN
    ZF(KSTART:KEND,0)       = PUF(KSTART:KEND,1)
    ZF(KSTART:KEND,KFLEV+1) = PUF(KSTART:KEND,KFLEV)
    ZF(KSTART:KEND,1:KFLEV) = PUF(KSTART:KEND,1:KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'FDER','00',KPROMA,KSTART,KEND,KFLEV,ZF,ZDUF)
    ZF(KSTART:KEND,0)       = PVF(KSTART:KEND,1)
    ZF(KSTART:KEND,KFLEV+1) = PVF(KSTART:KEND,KFLEV)
    ZF(KSTART:KEND,1:KFLEV) = PVF(KSTART:KEND,1:KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'FDER','00',KPROMA,KSTART,KEND,KFLEV,ZF,ZDVF)
  ELSE
    ZF(KSTART:KEND,0)       = 0.0_JPRB
    ZF(KSTART:KEND,KFLEV+1) = 0.0_JPRB
    ZF(KSTART:KEND,1:KFLEV) = PUF(KSTART:KEND,1:KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'FDER','11',KPROMA,KSTART,KEND,KFLEV,ZF,ZDUF)
    ZF(KSTART:KEND,1:KFLEV) = PVF(KSTART:KEND,1:KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'FDER','11',KPROMA,KSTART,KEND,KFLEV,ZF,ZDVF)
  ENDIF

  IF (PRESENT(PNHPPI)) THEN
    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PX(JROF,JLEV) =&
         & REAL(PNHPPI(JROF,JLEV),JPRD)/REAL(PRT(JROF,JLEV)*PLNPR(JROF,JLEV),JPRD)*&
         & REAL(PHIFL(JROF,JLEV)*ZDUF(JROF,JLEV)+PHIFM(JROF,JLEV)*ZDVF(JROF,JLEV),JPRD)&
         & /YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PX(JROF,JLEV) =&
         & 1._JPRD/REAL(PRT(JROF,JLEV)*PLNPR(JROF,JLEV),JPRD)*&
         & REAL(PHIFL(JROF,JLEV)*ZDUF(JROF,JLEV)+PHIFM(JROF,JLEV)*ZDVF(JROF,JLEV),JPRD)&
         & /YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      ENDDO
    ENDDO
  ENDIF

ELSE

  IF (PRESENT(PNHPPI)) THEN
    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PX(JROF,JLEV)=( &
         & REAL(PUH(JROF,JLEV)-PUF(JROF,JLEV),JPRD)*PHIHL(JROF,JLEV) + &
         & REAL(PVH(JROF,JLEV)-PVF(JROF,JLEV),JPRD)*PHIHM(JROF,JLEV) + &
         & REAL(PUF(JROF,JLEV)-PUH(JROF,JLEV-1),JPRD)*PHIHL(JROF,JLEV-1) + &
         & REAL(PVF(JROF,JLEV)-PVH(JROF,JLEV-1),JPRD)*PHIHM(JROF,JLEV-1) &
         & )*REAL(PNHPPI(JROF,JLEV),JPRD)/REAL(PRT(JROF,JLEV)*PLNPR(JROF,JLEV),JPRD)
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PX(JROF,JLEV)=( &
         & REAL(PUH(JROF,JLEV)-PUF(JROF,JLEV),JPRD)*PHIHL(JROF,JLEV) + &
         & REAL(PVH(JROF,JLEV)-PVF(JROF,JLEV),JPRD)*PHIHM(JROF,JLEV) + &
         & REAL(PUF(JROF,JLEV)-PUH(JROF,JLEV-1),JPRD)*PHIHL(JROF,JLEV-1) + &
         & REAL(PVF(JROF,JLEV)-PVH(JROF,JLEV-1),JPRD)*PHIHM(JROF,JLEV-1) &
         & )/REAL(PRT(JROF,JLEV)*PLNPR(JROF,JLEV),JPRD)
      ENDDO
    ENDDO
  ENDIF

ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPXX',1,ZHOOK_HANDLE)

END SUBROUTINE GPXX

