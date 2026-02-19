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

SUBROUTINE LAQMLIMITER(PRDETAR,LDQM3DCONS,KPROMA,KPROMB,KST,KPROF,KFLEV, &
 &                   KFLDN,KFLDX,KL0,PXSL,PXF)

! Purpose :
! -------
!   LAQMLIMITER - "aposteriori" quasi-monotone limiter for cubic-lagrange interpolation: 
!            Limits pre-computed PXF by surounding gp values. Final value never
!            exceeds min/max of surounding values and is contained in the interval
!            defined by the linear, cubic interpolated values. 
!            An option exists to distributed in excess/deficit in the vertical
!            improving thus mass conservation.

! Interface :
! ---------
!   INPUT:
!     PRDETAR - VRDETAR (vertical coordinate related par)
!     LDQM3DCONS - Bermejo & Staniforth quasi-monotone limiter, see type_gfld in yom_ygfl 
!     KPROMA  - horizontal dimension for grid-point quantities
!     KPROMB  - horizontal dimension for interpolation point quantities
!     KST     - first element of arrays where computations are performed
!     KPROF   - depth of work
!     KFLEV   - vertical dimension
!     KFLDN   - number of the first field
!     KFLDX   - number of the last field
!     KL0     - indices of the four western points of the 16 point
!               interpolation grid
!     PXSL    - quantity to be quasi-monotically limited 
!   OUTPUT:
!     PXF     - final limited interpolated variable

! Externals :
! ---------
!   None.

! Method :
! ------
!   Based on Bermejo, Staniforth MWR 1992, Vol. 20.

! Reference :
! ---------

! Author :
! ------
!   -Feb-2014 M. Diamantakis

! Modifications :
! -------------
! O. Marsden April 2017 : remove dependency on YOM_YGFL by passing in LDQM3DCONS as argument
! F. Vana   27-May-2020 : Simplified & use of precomputed (SP safe) RDETAR array

! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDETAR(KFLEV)
LOGICAL           ,INTENT(IN)    :: LDQM3DCONS
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KPROMA*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXF(KPROMB,KFLEV)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV1L1, IV1L2, IV2L1, IV2L2
INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZF111, ZF112, ZF121, ZF122
REAL(KIND=JPRB) :: ZF211, ZF212, ZF221, ZF222
REAL(KIND=JPRB) :: ZMIN, ZMAX, ZXF
REAL(KIND=JPRB) :: ZSURPL(KPROMB)
REAL(KIND=JPRB) :: ZRVETA(KFLEV)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LAQMLIMITER',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------

! offsets for getting values on interpolation stencil
!------------------------------------------------------------------------------
!  At LEV k=0,1,2,3, latitude line m:
!    KL0(,,m)+IVkL0: western-most from d.p. lev k point
!    KL0(,,m)+IVkL1: nearest western from d.p. NE lev k point
!    KL0(,,m)+IVkL2: nearest eastern from d.p. SW lev k point
!    KL0(,,m)+IVkL3: eastern-most from d.p. lev k point
!------------------------------------------------------------------------------
IV1L1=1+KPROMA
IV1L2=2+KPROMA
IV2L1=IV1L1+KPROMA
IV2L2=IV1L2+KPROMA
!------------------------------------------------------------------------------
! 8 point quasi-monotone correction (Bermejo & Staniforth)
!   with conservation (when activated)
!------------------------------------------------------------------------------

ZRVETA(1)=1.0_JPRB
ZRVETA(2:KFLEV)=PRDETAR(2:KFLEV)
ZSURPL(:)=0.0_JPRB

DO JLEV=1,KFLEV
  DO JROF=KST,KPROF
! compute surounding point values
    ZF111=PXSL(KL0(JROF,JLEV,1)+IV1L1) ! level above d.p. NW point
    ZF211=PXSL(KL0(JROF,JLEV,1)+IV1L2) ! level above d.p. NE point
    ZF121=PXSL(KL0(JROF,JLEV,2)+IV1L1) ! level above d.p. SW point
    ZF221=PXSL(KL0(JROF,JLEV,2)+IV1L2) ! level above d.p. SE point
    ZF112=PXSL(KL0(JROF,JLEV,1)+IV2L1) ! level below d.p. NW point
    ZF212=PXSL(KL0(JROF,JLEV,1)+IV2L2) ! level below d.p. NE point
    ZF122=PXSL(KL0(JROF,JLEV,2)+IV2L1) ! level below d.p. SW point
    ZF222=PXSL(KL0(JROF,JLEV,2)+IV2L2) ! level below d.p. SE point
! compute local min/max
    ZMIN=MIN(ZF111,ZF211,ZF121,ZF221,ZF112,ZF212,ZF122,ZF222)
    ZMAX=MAX(ZF111,ZF211,ZF121,ZF221,ZF112,ZF212,ZF122,ZF222)
! limit & store surplus
    ZXF=PXF(JROF,JLEV)
    IF (LDQM3DCONS ) ZXF=ZXF+ZSURPL(JROF)
    PXF(JROF,JLEV)=MAX(ZMIN,MIN(ZMAX,ZXF))
    IF (LDQM3DCONS ) ZSURPL(JROF)=(ZXF-PXF(JROF,JLEV))*ZRVETA(JLEV)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('LAQMLIMITER',1,ZHOOK_HANDLE)

END SUBROUTINE LAQMLIMITER

