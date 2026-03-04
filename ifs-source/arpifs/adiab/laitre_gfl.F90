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

SUBROUTINE LAITRE_GFL(YDGEOMETRY,YGFL,YDML_DYN,&
 & KGFL,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,KHVI,KWK,&
 & CDINT,PWENOALPHA,LDMOM,KNOWENO,KL0,PLSCAW,PRSCAW,PXSL,POUT,PXSPSL)

! Purpose :
! -------
!   LAITRE_GFL - semi-Lagrangian scheme: origin point interpolations 
!   for GFL variables.

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
!     KHVI       - 1/0: Cubic Hermite vert. interp. are needed/not needed
!     KWK        - set of horizontal non-linear weights to be used.
!     CDINT      - kind of interpolation used
!     PWENOALPHA - value of ALPHA (p) controlling level of non-oscillatory control
!     LDMOM      - momentum like variable versus heat (scalar) type of variable
!     KNOWENO    - special boundary treatment for WENO
!     KL0        - indices of the four western points of the 16 point
!                  interpolation grid
!     PLSCAW     - linear weights (distances) for interpolations.
!     PRSCAW     - non-linear weights for interpolations.
!     PXSL       - quantity to interpolate at O in U-wind eqn
!     PXSPSL     - spline representation of above
!   OUTPUT:
!     POUT       - interpolated quantity at O

! Externals :
! ---------
!   ABOR1
!   LAIHVT
!   LAITVSPCQM
!   LAITRI

! Method :
! ------
!   See documentation.

! Reference :
! ---------

! Author :
! ------
!   ??-Nov-2004, F. Vana (CHMI) and K. YESSAD (Meteo-France)

! Modifications :
! -------------
!   30-06-2008 J. Masek   Dataflow for new SLHD interpolators.
!   17-Sep-2008 F. Vana   Weights driven interp. for LAITVSPCQM and LAIHVT(..) 
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad Nov 2008: interpolation routines: merge QM with not-QM version.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   M. Diamantakis (June 2012): add code for quasi-monotone mass fixer 
!   M. Diamantakis (Dec  2012): add code for additional mass fixers
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   F. Vana   13-Feb-2014  Distinguish SLHD for momentum and heat variables.
!   M. Diamantakis (Feb 2014): - introduce new QM interpol limiter
!                              - make code independent on use of mass fixer type
!   M. Diamantakis (Feb 2016): - introduce new LQML3D interpol limiter
!                              - improve efficiency of LQM3D
!   F. Vana  20-Feb-2019  Quintic vertical interpolation
!   M. Diamantakis & F. Vana May-2020: LQM3D for cubic ans quintic interpolation
! End Modifications
!   -----------------------------------------------------------------------

USE YOM_YGFL     , ONLY : TYPE_GFLD
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE

!   -----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KHVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWK
CHARACTER(LEN=12) ,INTENT(IN)    :: CDINT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWENOALPHA
LOGICAL           ,INTENT(IN)    :: LDMOM
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(KPROMA,KFLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMA,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAW(KPROMA,KFLEV,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAW(KPROMA,KFLEV,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL    :: PXSPSL(KASLB1,KFLDN:KFLDX) 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!   -----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "laihvt.intfb.h"
#include "laitvspcqm.intfb.h"
#include "laitri.intfb.h"
#include "laitriqm3d.intfb.h"
#include "laqmlimiter.intfb.h"
#include "laitri_weno.intfb.h"

!   -----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAITRE_GFL',0,ZHOOK_HANDLE)

ASSOCIATE(YDTLSCAW=>YDML_DYN%YYTLSCAW,YDTRSCAW=>YDML_DYN%YYTRSCAW,YDDYN=>YDML_DYN%YRDYN, & 
  &       YDVSLETA=>YDML_DYN%YRSLINT%YRVSLETA,YDVETA=>YDGEOMETRY%YRVETA, LQM3DCONS=>YGFL%LQM3DCONS)
!   -----------------------------------------------------------------------

! 1. Origin point interpolations for GFL variables
!    ---------------------------------------------

! 1.1 SLHD interpolations
!     -------------------
! NOTE:
!   The SLHD versions of the LAIHVT and
!   LAITVSPCQM interpolations are currently not coded.

IF (CDINT == 'LAITQM(SLD) ') THEN

  IF (LDMOM) THEN
    ! Lagrange cubic polynomial with SLHD (momentum), quasi-monotonic
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLD(KWK)),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLD(KWK)),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLD),&
     & PXSL,POUT)
  ELSE
    ! Lagrange cubic polynomial with SLHD (heat), quasi-monotonic
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLT),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLT),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLT),&
     & PXSL,POUT)
  ENDIF

ELSEIF (CDINT == 'LAITQMH(SLD)') THEN

  IF (LDMOM) THEN
    ! Lagrange cubic polynomial with SLHD (momentum), quasi-monotonic in horizontal
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLD(KWK)),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLD(KWK)),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLD),&
     & PXSL,POUT)
  ELSE
    ! Lagrange cubic polynomial with SLHD (heat), quasi-monotonic in horizontal
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLT),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLT),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLT),&
     & PXSL,POUT)
  ENDIF

ELSEIF (CDINT == 'LAITRI(SLD) ') THEN

  IF (LDMOM) THEN
    ! Lagrange cubic polynomial with SLHD (momentum)
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLD(KWK)),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLD(KWK)),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLD),&
     & PXSL,POUT)
  ELSE
    ! Lagrange cubic polynomial with SLHD (heat)
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLT),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLT),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLT),&
     & PXSL,POUT)
  ENDIF

ELSEIF (CDINT == 'LAITQM(MAD) ') THEN

! 1.2 COMAD (cubic) interpolations

  ! Lagrange cubic polynomial with COMAD, quasi-monotonic
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PRSCAW(1,1,YDTRSCAW%M_WCLAMAD(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLOMAD),PRSCAW(1,1,YDTRSCAW%M_WCLOMAD(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWMAD),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAITQM3(MAD)') THEN

  ! Lagrange cubic polynomial, quasi-monotonic, 3D Bermejo & Staniforth limiter
  CALL LAITRIQM3D(YDGEOMETRY%YRVETA,YGFL%LQM3DCONS,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
    & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PRSCAW(1,1,YDTRSCAW%M_WCLAMAD(KWK)),&
    & PLSCAW(1,1,YDTLSCAW%M_WDLOMAD),PRSCAW(1,1,YDTRSCAW%M_WCLOMAD(KWK)),&
    & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWMAD),PLSCAW(1,1,YDTLSCAW%M_WDVERMAD),&
    & PXSL,POUT)

ELSEIF (CDINT == 'LAITQML(MAD)') THEN

  ! Lagrange cubic polynomial, quasi-monotonic, 3D Bermejo & Staniforth limiter
  CALL LAITRIQM3D(YDGEOMETRY%YRVETA,YGFL%LQM3DCONS,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
    & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PRSCAW(1,1,YDTRSCAW%M_WCLAMAD(KWK)),&
    & PLSCAW(1,1,YDTLSCAW%M_WDLOMAD),PRSCAW(1,1,YDTRSCAW%M_WCLOMAD(KWK)),&
    & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWMAD),PLSCAW(1,1,YDTLSCAW%M_WDVERMAD),&
    & PXSL,POUT)

ELSEIF (CDINT == 'LAITQM3D(SLD') THEN
  ! NOTE: If one gets here, then the abort line can be commented out. 
  !       For the moment the GFL setup doesn't know this possibility.
  CALL ABOR1('LAITRE_GFL: LAITQM3D(SLD IS NOT YET AVAILABLE')

  IF (LDMOM) THEN
    ! Lagrange cubic polynomial with SLHD (momentum)
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLD(KWK)),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLD(KWK)),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLD),&
     & PXSL,POUT)
  ELSE
    ! Lagrange cubic polynomial with SLHD (heat)
    CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
     & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLASLT),&
     & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLOSLT),&
     & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWSLT),&
     & PXSL,POUT)
  ENDIF
  CALL LAQMLIMITER(YDVSLETA%VRDETAR,LQM3DCONS,&
   & KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,KL0,PXSL,POUT)


ELSEIF (CDINT == 'LAITQMH(MAD)') THEN
  ! Lagrange cubic polynomial with COMAD, quasi-monotonic in horizontal
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PRSCAW(1,1,YDTRSCAW%M_WCLAMAD(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLOMAD),PRSCAW(1,1,YDTRSCAW%M_WCLOMAD(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWMAD),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAITRI(MAD) ') THEN

  ! Lagrange cubic polynomial with COMAD
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAMAD),PRSCAW(1,1,YDTRSCAW%M_WCLAMAD(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLOMAD),PRSCAW(1,1,YDTRSCAW%M_WCLOMAD(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWMAD),&
   & PXSL,POUT)

! 1.3 non-SLHD non-COMAD interpolations
!     -----------------------

ELSEIF (CDINT == 'LAIHVTQM    ') THEN 

  ! Hermite cubic polynomial on the vertical "HVT", quasi-monotonic
  CALL LAIHVT(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVDERW),PRSCAW(1,1,YDTRSCAW%M_WHVW),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAIHVTQMH   ') THEN

  ! Hermite cubic polynomial on the vertical "HVT", quasi-monotonic 
  ! in horizontal
  CALL LAIHVT(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVDERW),PRSCAW(1,1,YDTRSCAW%M_WHVW),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAIHVT      ') THEN

  ! Hermite cubic polynomial on the vertical "HVT"
  CALL LAIHVT(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVDERW),PRSCAW(1,1,YDTRSCAW%M_WHVW),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAITVSPCQM  ' .AND. PRESENT(PXSPSL)) THEN

  ! spline "VSPC", quasi-monotonic
  CALL LAITVSPCQM(YDVSLETA,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTWS),&
   & PXSPSL,POUT,PXSL)

ELSEIF (CDINT == 'LAITVSPCQM  ' .AND. .NOT.PRESENT(PXSPSL)) THEN

        CALL ABOR1('LAITRE_GFL: MISSING ARGUMENT PXSPSL FOR INTERPOLATION TYPE : ' // CDINT)

ELSEIF (CDINT == 'LAITQM3D    ') THEN

  ! Lagrange cubic polynomial, quasi-monotonic, 3D Bermejo & Staniforth limiter
  CALL LAITRIQM3D(YDVETA,LQM3DCONS,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
    & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
    & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
    & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTW),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
    & PXSL,POUT)

 !! Alternative bit-reproducible method consistent with TL/AD code
 !CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,-1,&
 ! & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
 ! & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
 ! & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTW),&
 ! & PXSL,POUT,PRDETAR=YDVSLETA%VRDETAR,LDQM3DCONS=LQM3DCONS)

  ! Previous could be written in a more general way using two step method
  !  calling 1/ any interpolation routine and 2/ aposteriori LAQMLIMITER fixer.
  !  This alternative is not very practical for adjoint code as it implies
  !  a three step sequence there: 1/Trajectory interpolation 2/ LAQMLIMITER fixer
  !  3/ Adjoint interpolation (with partial trajectory again)

ELSEIF (CDINT == 'LAITQML3D   ') THEN

  ! Lagrange cubic polynomial, quasi-monotonic, 3D linear interp limiter
  CALL LAITRIQM3D(YDVETA,LQM3DCONS,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
    & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
    & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
    & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTW),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
    & PXSL,POUT)

ELSEIF (CDINT == 'LAITQM      ') THEN

  ! Lagrange cubic polynomial, quasi-monotonic
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
    & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
    & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
    & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTW),&
    & PXSL,POUT)

ELSEIF (CDINT == 'LAITQMH     ') THEN

  ! Lagrange cubic polynomial, quasi-monotonic in horizontal
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTW),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAITRI      ') THEN

  ! Lagrange cubic polynomial
  CALL LAITRI(KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)),&
   & PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)),&
   & KL0,PRSCAW(1,1,YDTRSCAW%M_WVINTW),&
   & PXSL,POUT)

ELSEIF (CDINT == 'LAITRWENOQM ') THEN

  ! WENO quintic on the vertical, quasi-monotonic
  CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,2,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)), &
   & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)), &
   & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,YDTRSCAW%M_WVINTW), &
   & PXSL,POUT,PALPHA=PWENOALPHA)

ELSEIF (CDINT == 'LAITRWENOQM3') THEN

  ! WENO quintic on the vertical with Bermejo & Staniforth fixer
  CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,-1,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)), &
   & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)), &
   & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,YDTRSCAW%M_WVINTW), &
   & PXSL,POUT,PALPHA=PWENOALPHA,PRDETAR=YDVSLETA%VRDETAR,LDQM3DCONS=LQM3DCONS)

ELSEIF (CDINT == 'LAITRWENOQMH') THEN

  ! WENO quintic on the vertical, horizontally quasi-monotonic
  CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,1,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)), &
   & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)), &
   & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,YDTRSCAW%M_WVINTW), &
   & PXSL,POUT,PALPHA=PWENOALPHA)

ELSEIF (CDINT == 'LAITRWENO   ') THEN

  ! WENO quintic on the vertical
  CALL LAITRI_WENO(YDDYN,KASLB1,KPROMA,KSTART,KPROF,KFLEV,KFLDN,KFLDX,0,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PRSCAW(1,1,YDTRSCAW%M_WCLA(KWK)), &
   & PLSCAW(1,1,YDTLSCAW%M_WDLO) ,PRSCAW(1,1,YDTRSCAW%M_WCLO(KWK)), &
   & KL0,KNOWENO,PRSCAW(1,1,YDTRSCAW%M_CW),PRSCAW(1,1,YDTRSCAW%M_WVINTW), &
   & PXSL,POUT,PALPHA=PWENOALPHA)

ELSE
  CALL ABOR1('LAITRE_GFL: UNKNOWN INTERPOLATION TYPE: ' // CDINT)
ENDIF

!   -----------------------------------------------------------------------
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('LAITRE_GFL',1,ZHOOK_HANDLE)
END SUBROUTINE LAITRE_GFL
