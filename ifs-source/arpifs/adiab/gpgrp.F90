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
SUBROUTINE GPGRP(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,KST,KEND,&
 & PRT,PRTL,PRTM,PREL,PREM,PXYB,PXYBDER,&
 & PHIHL,PHIHM,PHIFL,PHIFM,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PSGRTL,PSGRTM)

!**** *GPGRP* - Computation of the pressure gradient force term used in the
!               RHS of the horizontal wind equation in hydrostatic model.
!               This term is computed at full levels.

!     Purpose.
!     --------

!        The pressure gradient force term writes:
!          - grad[gz] - RT grad[log(prehyds)]

!**   Interface.
!     ----------
!        *CALL* *GPGRP(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KST         : start of work
!           KEND        : working length
!           PRT         : "R(air)*temperature" at full levels
!           PRTL        : zonal component of "grad (RT)" at full levels
!           PRTM        : meridian component of "grad (RT)" at full levels
!           PREL        : zonal component of "grad (prehyds)"
!           PREM        : meridian component of "grad (prehyds)"
!           PXYB        : contains pressure depth, "delta", "alpha".
!           PXYBDER     : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
!           PHIHL       : zonal component of "grad[gz]" at half levels
!           PHIHM       : merid component of "grad[gz]" at half levels
!           PHIFL       : zonal component of "grad[gz]" at full levels
!           PHIFM       : merid component of "grad[gz]" at full levels

!         * OUTPUT:
!           PSGRTL      : zonal comp. of pressure gradient force
!           PSGRTM      : merid comp. of pressure gradient force

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See  documentation Arpege ALGORITHME CHAPITRE 6 paragraphe 6

!     Author.
!     -------
!      Florence Rabier
!      Original : 91-11-01

!     Modifications.
!     --------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   N. Wedi (Nov 2011) : cleaning
!   K. Yessad (Mar 2017): move NHEE code in new GNHEE_GRP.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE INTDYN_MOD   , ONLY : YYTXYB, YYTXYBDER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYBDER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIFM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!!! #include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,NPROMA=>YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF THE PRESSURE GRADIENT TERM.
!              ------------------------------------------

IF (YDGEOMETRY%YRCVER%LVERTFE) THEN

  ! * The pressure gradient term is directly computed at full levels in
  !   this case and no quantity at half layers is computed.

  ! * First transfer " grad[gz] + RT [grad (log(prehyd))] " in (PSGRTL,PSGRTM).
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PHIFL(JROF,JLEV) &
       & +PRT(JROF,JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RTGR)*PREL(JROF)
      PSGRTM(JROF,JLEV)=PHIFM(JROF,JLEV) &
       & +PRT(JROF,JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RTGR)*PREM(JROF)
    ENDDO
  ENDDO

ELSE

  ! * The pressure gradient term calculation requires some
  !   calculations on half levels.

  ! * In the hydrostatic model, an accurate way to code this
  !   is to exhibit the quantity "grad alpha + log prehyd"
  !   which is simpler to compute than the separate terms
  !   "grad alpha" and "log prehyd"
  !   (some quantities cancel each other).
  !   One starts from half level "grad(Phi)" and write:
  !   [grad(Phi) + RT grad(log(prehyd))]_l
  !   = [grad(Phi)]_lbar
  !   + [alpha]_l [grad(RT)]_l
  !   + [RT]_l [grad(alpha + log(prehyd))]_l

  ! * First get "grad[gz]" on half levels.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PHIHL(JROF,JLEV)
      PSGRTM(JROF,JLEV)=PHIHM(JROF,JLEV)
    ENDDO
  ENDDO

  ! * Calculation of [ grad(Phi) + RT grad(log(prehyd)) ]
  !   at full levels (cf. what is done in GPGRGEO to compute "grad(Phi)"
  !   at full levels, but replacing (alphl,alphm) by (alphpll,alphplm)):

  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PSGRTL(JROF,JLEV)&
       & +(PXYB(JROF,JLEV,YYTXYB%M_ALPH)*PRTL(JROF,JLEV)&
       & +PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)*PRT(JROF,JLEV))
      PSGRTM(JROF,JLEV)=PSGRTM(JROF,JLEV)&
       & +(PXYB(JROF,JLEV,YYTXYB%M_ALPH)*PRTM(JROF,JLEV)&
       & +PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)*PRT(JROF,JLEV))
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPGRP',1,ZHOOK_HANDLE)
END SUBROUTINE GPGRP
