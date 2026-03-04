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

SUBROUTINE LAITRE_GMV( &
 ! --- INPUT ---------------------------------------------------------
 & YDML_DYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,KWK, &
 & LDLIN,LDCOMAD,LDSLHD,LDVWENO,PWENOALPHA,LDMOM,LDQM,LDQMH,KNOWENO,&
 & KL0,PLSCAW,PRSCAW,PXSL, &
 ! --- OUTPUT --------------------------------------------------------
 & POUT)

! Purpose :
! -------
!   LAITRE_GMV - semi-Lagrangian scheme: origin point interpolations 
!   for GMV and some miscellaneous variables.

! Interface :
! ---------
!   INPUT:
!     KASLB1     - horizontal dimension of SLBUF1
!     KPROMA     - horizontal dimension
!     KSTART     - first element of arrays where computations are performed
!     KPROF      - depth of work
!     KFLEV      - vertical dimension
!     KFLDN      - number of the first field
!     KFLDX      - number of the last field
!     KWK        - set of horizontal non-linear weights to be used.
!     LDLIN      - switch for tri-linear interpolation
!     LDCOMAD    - switch for COMAD interpolation
!     LDSLHD     - switch: .T. = diffusive interpolation (SLHD)
!                          .F. = normal interpolation
!     LDVWENO    - if .T. vertical WENO is used with PWENOALPHA
!     PWENOALPHA - value of ALPHA (p in literature) defining level
!                  of non-oscillatory control
!     LDMOM      - if .T. momentum variable, otherwise heat variable
!     LDQM       - use three-dimensional quasi-monotone interpolation
!     LDQMH      - use quasi-monotone interpolation in the horizontal
!     KNOWENO    - special boundary treatment for WENO
!     KL0        - indices of the four western points of the 16 point 
!                  interpolation grid
!     PLSCAW     - linear weights (distances) for interpolations.
!     PRSCAW     - non-linear weights for interpolations.
!     PXSL       - quantity to interpolate at O
!   OUTPUT:
!     POUT       - interpolated quantity at O

! Externals :
! ---------
!   LAITRI

! Method :
! ------
!   See documentation.

! Reference :
! ---------

! Author :
! ------
!   ??-Nov-2004 F. Vana (CHMI) and K. Yessad (Meteo-France)

! Modifications :
! -------------
!    28-Aug-2007  F. Vana  cleaning out the 4-points splines
!    30-06-2008 J. Masek   Dataflow for new SLHD interpolators.
!    K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!    K. Yessad Nov 2008: interpolation routines: merge QM with not-QM version.
!    K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!    M. Diamantakis (May 2012): add option for Bermejo & Staniforth QMSL filter
!                               for consistency with the quasi-monotone mass fixer 
!    S. Malardel (Nov 2013): COMAD weights for SL interpolations
!    F. Vana 13-Feb-2014 Distinguish between momentum and heat SLHD.
!    K. Yessad (June 2017): introduce LDLIN.
!    F. Vana  20-Feb-2019 Quintic vertical interpolation 
!    F. Vana  18-Jul-2019:  SLVF
! End Modifications
!    ----------------------------------------------------------------------

USE PARKIND1    , ONLY : JPIM     ,JPRB
USE YOMHOOK     , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE

!    ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWK
LOGICAL           ,INTENT(IN)    :: LDLIN
LOGICAL           ,INTENT(IN)    :: LDCOMAD
LOGICAL           ,INTENT(IN)    :: LDSLHD
LOGICAL           ,INTENT(IN)    :: LDVWENO
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWENOALPHA
LOGICAL           ,INTENT(IN)    :: LDMOM
LOGICAL           ,INTENT(IN)    :: LDQM
LOGICAL           ,INTENT(IN)    :: LDQMH
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(KPROMA,KFLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMA,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAW(KPROMA,KFLEV,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAW(KPROMA,KFLEV,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(KPROMA,KFLEV) 

!    ----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IQM, IM_VERT
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!    ----------------------------------------------------------------------

#include "laitli.intfb.h"
#include "laitri.intfb.h"
#include "laitri_weno.intfb.h"

!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAITRE_GMV',0,ZHOOK_HANDLE)

ASSOCIATE(YDTLSCAW=>YDML_DYN%YYTLSCAW,YDTRSCAW=>YDML_DYN%YYTRSCAW,YDDYN=>YDML_DYN%YRDYN, &
 & YDDYNA=>YDML_DYN%YRDYNA)
!    ----------------------------------------------------------------------

! 0. Setup
!    -----

! LSLVF (only applies to temperature and WENO or Lagrangian cubic interpolation)
IF ((.NOT.LDMOM).AND.YDDYNA%LSLVF) THEN
  IM_VERT=YDTRSCAW%M_WVINTWSLVF
ELSE
  IM_VERT=YDTRSCAW%M_WVINTW
ENDIF

IF (LDQM) THEN
  IQM=2
ELSEIF (LDQMH) THEN
  IQM=1
ELSE
  IQM=0
ENDIF

! 1. Origin point interpolations for GMV and some miscellaneous variables
!    --------------------------------------------------------------------

IF (LDLIN) THEN
  IF (LDCOMAD) THEN
    CALL LAITLI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PLSCAW(1,1,YDTLSCAW%M_WDLOMAD+1),&
     & KL0(1,1,1),PLSCAW(1,1,YDTLSCAW%M_WDVERMAD),&
     & PXSL,POUT)
  ELSE
    CALL LAITLI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PLSCAW(1,1,YDTLSCAW%M_WDLO+1),&
     & KL0(1,1,1),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
     & PXSL,POUT)
  ENDIF
ELSEIF (LDSLHD) THEN
  ! NOTE: When both WENO and SLHD are set WENO takes over for the vertical part.
  !       This tricky part is usefull for the upper 20hPa treatment.
  IF (LDMOM) THEN
    ! Momentum variables
    IF (LDVWENO) THEN
      CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
       & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLD(KWK)), &
       & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLOSLD(KWK)), &
       & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,YDTRSCAW%M_WVINTW), &
       & PXSL,POUT,PALPHA=PWENOALPHA)
    ELSE
      CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
       & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLD(KWK)),&
       & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLD(KWK)),&
       & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLD),&
       & PXSL,POUT)
    ENDIF
  ELSE
    ! Heat variables
    IF (LDVWENO) THEN
      CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
       & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLT), &
       & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLOSLT), &
       & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,IM_VERT), &
       & PXSL,POUT,PALPHA=PWENOALPHA)
    ELSE
      CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
       & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLT),&
       & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLT),&
       & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLT),&
       & PXSL,POUT)
    ENDIF
  ENDIF
ELSEIF (LDCOMAD) THEN
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PRSCAW(1,1,YDTRSCAW%M_WCLAMAD(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLOMAD),PRSCAW(1,1,YDTRSCAW%M_WCLOMAD(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWMAD),&
   & PXSL,POUT)
ELSEIF (LDVWENO) THEN
  CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)), &
   & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)), &
   & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,IM_VERT), &
   & PXSL,POUT,PALPHA=PWENOALPHA)
ELSE
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,IQM,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,IM_VERT),&
   & PXSL,POUT)
ENDIF

!    ----------------------------------------------------------------------
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('LAITRE_GMV',1,ZHOOK_HANDLE)
END SUBROUTINE LAITRE_GMV
