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

!option! -O extendreorder
SUBROUTINE LAITRIQM3D(YDVETA,LDQM3DCONS,KSLB1,KPROMA,KST,KPROF,KFLEV, &
 & KFLDN,KFLDX,KQM,PDLAT,PCLA,PDLO,PCLO,KL0,PVINTW, &
 & PDVER,PXSL,PXF)

! Purpose :
! -------
!   LAITRIQM3D - semi-Lagrangian scheme: tri-dimensional 32-point
!   interpolations with optional 3d quasi-monotonic limiters built in.
!   KQM=1 is a combined ecmwf quasi-cubic/linear scheme which is by 
!   construction quasi-monotone while KQM=2 is the standard 
!   Bermejo & Staniforth applied in 3D. 
!   KQM=2 is equivavalent to calling LAITRI(KQM=0) with
!   and then LAQMLITER() but slightly more efficient here.
!   Differs form the standard limiter KQM=2 in laitri as it is 3d
!   limiter applied only once at the end of all 1-d interpolations. 
!    

! Interface :
! ---------
!   INPUT:
!     YDVETA - defining the vertical coordinate: eta
!     LDQM3DCONS - Bermejo & Staniforth quasi-monotone limiter, see type_gfld in yom_ygfl 
!     KSLB1   - horizontal dimension for grid-point quantities
!     KPROMA  - horizontal dimension for interpolation point quantities
!     KST     - first element of arrays where computations are performed
!     KPROF   - depth of work
!     KFLEV   - vertical dimension
!     KFLDN   - number of the first field
!     KFLDX   - number of the last field
!     KQM     - type of limiter
!               1: standard Bermejo-Staniforth
!               2: use linear interpolant which is qm by definition
!     PDLAT   - distance for horizontal linear interpolations in latitude
!     PCLA    - weights for horizontal cubic interpolations in latitude
!     PDLO    - distances for horizontal linear interpolations
!               in longitude (latitude rows 0, 1, 2, 3)
!     PCLO    - weights for horizontal cubic interpolations in longitude 
!               (latitude rows 1, 2)
!     KL0     - indices of the four western points of the 16 point
!               interpolation grid
!     PVINTW  - weights for cubic vertical interpolation
!     PDVER   - weights (distances) for vertical linear interpolation
!                    on a same vertical.
!     PXSL    - quantity to be interpolated
!   OUTPUT:
!     PXF     - interpolated variable

! Externals :
! ---------
!   None.

! Method :
! ------
!   cubic interpolation with embedded 3D quasi-monotone limiter. Two oprions are available:
!     - LQM3D limiter contained in LAQMLIMITER() but coded in a more efficient way
!     - LQML3D limiter i.e. use the linear interpolant when a new min/max generated    

! Reference :
! ---------

! Author :
! ------
! Michail Diamantakis

! Modifications :
! -------------
! Original : Feb-2016
! O. Marsden April 2017 : remove dependency on YOM_YGFL by passing in LDQM3DCONS as argument
! End modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPIA
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVERT  , ONLY : TVETA

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVETA)       ,INTENT(IN)  :: YDVETA
LOGICAL           ,INTENT(IN)  :: LDQM3DCONS
INTEGER(KIND=JPIM),INTENT(IN)  :: KSLB1 
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)  :: KST 
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)  :: KQM
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDLAT(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCLA(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDLO(KPROMA,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCLO(KPROMA,KFLEV,3,2)
INTEGER(KIND=JPIM),INTENT(IN)  :: KL0(KPROMA,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PVINTW(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDVER(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PXSL(KSLB1*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(OUT) :: PXF(KPROMA,KFLEV)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV0L1, IV0L2, IV1L0, IV1L1, IV1L2, IV1L3, &
 & IV2L0, IV2L1, IV2L2, IV2L3, IV3L1, IV3L2
INTEGER(KIND=JPIM) :: JLEV, JROF, IQM

REAL(KIND=JPRB) :: ZF111(KPROMA), ZF112(KPROMA), ZF121(KPROMA), ZF122(KPROMA)
REAL(KIND=JPRB) :: ZF211(KPROMA), ZF212(KPROMA), ZF221(KPROMA), ZF222(KPROMA)
REAL(KIND=JPRB) :: Z01(KPROMA), Z02(KPROMA)
REAL(KIND=JPRB) :: Z10(KPROMA), Z11(KPROMA), Z12(KPROMA), Z13(KPROMA)
REAL(KIND=JPRB) :: Z20(KPROMA), Z21(KPROMA), Z22(KPROMA), Z23(KPROMA)
REAL(KIND=JPRB) :: Z31(KPROMA), Z32(KPROMA)
REAL(KIND=JPRB) :: ZRVETA(KFLEV)
REAL(KIND=JPRB) :: ZSURPL(KPROMA)
REAL(KIND=JPRB) :: Z0, Z1, Z2, Z3, ZMIN, ZMAX, ZDMAX, ZDMIN, ZXF
REAL(KIND=JPRB) :: ZINF, ZINFLO1, ZINFLO2, ZSUP, ZSUPLO1, ZSUPLO2 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAITRIQM3D',0,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

! 1. Interpolations
!    --------------

! offsets for getting values on interpolation stencil
IV0L1=1
IV0L2=2
IV1L0=  KSLB1
IV1L1=1+KSLB1
IV1L2=2+KSLB1
IV1L3=3+KSLB1
IV2L0=IV1L0+KSLB1
IV2L1=IV1L1+KSLB1
IV2L2=IV1L2+KSLB1
IV2L3=IV1L3+KSLB1
IV3L1=IV2L1+KSLB1
IV3L2=IV2L2+KSLB1

! 8 point quasi-monotone correction with vertical adjustment
! to improve conservation
IF (LDQM3DCONS.AND.(KQM==2)) THEN
  IQM=3
  ZRVETA(1)=1.0_JPRB
  DO JLEV=2,KFLEV
    ZRVETA(JLEV)=(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))/(YDVETA%VETAH(JLEV-1)-YDVETA%VETAH(JLEV-2))
  ENDDO
  ZSURPL(:)=0.0_JPRB
ELSE
  IQM=KQM
ENDIF

! 32-point interpolations
DO JLEV=1,KFLEV
    ! interpolations in longitude
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    ! interpolations in longitude, stencil level 0
    Z10(JROF)=PXSL(KL0(JROF,JLEV,1)+IV0L1)+PDLO(JROF,JLEV,1)&
      & *(PXSL(KL0(JROF,JLEV,1)+IV0L2)-PXSL(KL0(JROF,JLEV,1)+IV0L1))
    Z20(JROF)=PXSL(KL0(JROF,JLEV,2)+IV0L1)+PDLO(JROF,JLEV,2)&
      & *(PXSL(KL0(JROF,JLEV,2)+IV0L2)-PXSL(KL0(JROF,JLEV,2)+IV0L1))
    ! interpolations in longitude, stencil level 1
    Z01(JROF)=PXSL(KL0(JROF,JLEV,0)+IV1L1)+PDLO(JROF,JLEV,0)&
      & *(PXSL(KL0(JROF,JLEV,0)+IV1L2)-PXSL(KL0(JROF,JLEV,0)+IV1L1))
  ENDDO
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    ZF111(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L1)
    ZF211(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L2)
    Z11(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L0)&
      & +PCLO(JROF,JLEV,1,1)*(ZF111(JROF)-PXSL(KL0(JROF,JLEV,1)+IV1L0))&
      & +PCLO(JROF,JLEV,2,1)*(ZF211(JROF)-PXSL(KL0(JROF,JLEV,1)+IV1L0))&
      & +PCLO(JROF,JLEV,3,1)*&
      & (PXSL(KL0(JROF,JLEV,1)+IV1L3)-PXSL(KL0(JROF,JLEV,1)+IV1L0))
  ENDDO
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    ZF121(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L1)
    ZF221(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L2)
    Z21(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L0)&
      & +PCLO(JROF,JLEV,1,2)*(ZF121(JROF)-PXSL(KL0(JROF,JLEV,2)+IV1L0))&
      & +PCLO(JROF,JLEV,2,2)*(ZF221(JROF)-PXSL(KL0(JROF,JLEV,2)+IV1L0))&
      & +PCLO(JROF,JLEV,3,2)*&
      & (PXSL(KL0(JROF,JLEV,2)+IV1L3)-PXSL(KL0(JROF,JLEV,2)+IV1L0))
  ENDDO
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    Z31(JROF)=PXSL(KL0(JROF,JLEV,3)+IV1L1)+PDLO(JROF,JLEV,3)&
      & *(PXSL(KL0(JROF,JLEV,3)+IV1L2)-PXSL(KL0(JROF,JLEV,3)+IV1L1))
    ! interpolations in longitude, stencil level 2
    Z02(JROF)=PXSL(KL0(JROF,JLEV,0)+IV2L1)+PDLO(JROF,JLEV,0)&
      & *(PXSL(KL0(JROF,JLEV,0)+IV2L2)-PXSL(KL0(JROF,JLEV,0)+IV2L1))
  ENDDO
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    ZF112(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L1)
    ZF212(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L2)
    Z12(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L0)&
      & +PCLO(JROF,JLEV,1,1)*(ZF112(JROF)-PXSL(KL0(JROF,JLEV,1)+IV2L0))&
      & +PCLO(JROF,JLEV,2,1)*(ZF212(JROF)-PXSL(KL0(JROF,JLEV,1)+IV2L0))&
      & +PCLO(JROF,JLEV,3,1)*&
      & (PXSL(KL0(JROF,JLEV,1)+IV2L3)-PXSL(KL0(JROF,JLEV,1)+IV2L0))
  ENDDO
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    ZF122(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L1)
    ZF222(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L2)
    Z22(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L0)&
      & +PCLO(JROF,JLEV,1,2)*(ZF122(JROF)-PXSL(KL0(JROF,JLEV,2)+IV2L0))&
      & +PCLO(JROF,JLEV,2,2)*(ZF222(JROF)-PXSL(KL0(JROF,JLEV,2)+IV2L0))&
      & +PCLO(JROF,JLEV,3,2)*&
      & (PXSL(KL0(JROF,JLEV,2)+IV2L3)-PXSL(KL0(JROF,JLEV,2)+IV2L0))
  ENDDO
!DIR$ PREFERVECTOR
  DO JROF=KST,KPROF
    Z32(JROF)=PXSL(KL0(JROF,JLEV,3)+IV2L1)+PDLO(JROF,JLEV,3)&
      & *(PXSL(KL0(JROF,JLEV,3)+IV2L2)-PXSL(KL0(JROF,JLEV,3)+IV2L1))
    ! interpolations in longitude, stencil level 3
    Z13(JROF)=PXSL(KL0(JROF,JLEV,1)+IV3L1)+PDLO(JROF,JLEV,1)&
      & *(PXSL(KL0(JROF,JLEV,1)+IV3L2)-PXSL(KL0(JROF,JLEV,1)+IV3L1))
    Z23(JROF)=PXSL(KL0(JROF,JLEV,2)+IV3L1)+PDLO(JROF,JLEV,2)&
      & *(PXSL(KL0(JROF,JLEV,2)+IV3L2)-PXSL(KL0(JROF,JLEV,2)+IV3L1))
  ENDDO

  DO JROF=KST,KPROF
    ! interpolations in latitude, stencil levels 0, 1, 2, 3
    Z0=Z10(JROF)+PDLAT(JROF,JLEV)*(Z20(JROF)-Z10(JROF))
    Z1=Z01(JROF)+PCLA(JROF,JLEV,1)*(Z11(JROF)-Z01(JROF))&
      & +PCLA(JROF,JLEV,2)*(Z21(JROF)-Z01(JROF))&
      & +PCLA(JROF,JLEV,3)*(Z31(JROF)-Z01(JROF))
    Z2=Z02(JROF)+PCLA(JROF,JLEV,1)*(Z12(JROF)-Z02(JROF))&
      & +PCLA(JROF,JLEV,2)*(Z22(JROF)-Z02(JROF))&
      & +PCLA(JROF,JLEV,3)*(Z32(JROF)-Z02(JROF))
    Z3=Z13(JROF)+PDLAT(JROF,JLEV)*(Z23(JROF)-Z13(JROF))
    ! final interpolation in vertical
    PXF(JROF,JLEV)=Z0&
      & +PVINTW(JROF,JLEV,1)*(Z1-Z0)&
      & +PVINTW(JROF,JLEV,2)*(Z2-Z0)&
      & +PVINTW(JROF,JLEV,3)*(Z3-Z0)  
  ENDDO
  
  IF (IQM==1) THEN
    DO JROF=KST,KPROF
      ! compute local min/max
      ZMIN=MIN(ZF111(JROF),ZF211(JROF),ZF121(JROF),ZF221(JROF),&
        &       ZF112(JROF),ZF212(JROF),ZF122(JROF),ZF222(JROF))
      ZMAX=MAX(ZF111(JROF),ZF211(JROF),ZF121(JROF),ZF221(JROF),&
        &       ZF112(JROF),ZF212(JROF),ZF122(JROF),ZF222(JROF))
      ! limit
      IF (PXF(JROF,JLEV)>ZMAX.OR.PXF(JROF,JLEV)<ZMIN) THEN
        !     Linear interpolation.
        ZSUPLO1 = ZF111(JROF)+PDLO(JROF,JLEV,1)*(ZF211(JROF)-ZF111(JROF))
        ZSUPLO2 = ZF121(JROF)+PDLO(JROF,JLEV,2)*(ZF221(JROF)-ZF121(JROF))  
        ZINFLO1 = ZF112(JROF)+PDLO(JROF,JLEV,1)*(ZF212(JROF)-ZF112(JROF))  
        ZINFLO2 = ZF122(JROF)+PDLO(JROF,JLEV,2)*(ZF222(JROF)-ZF122(JROF))
        ZSUP    = ZSUPLO1+PDLAT(JROF,JLEV)*(ZSUPLO2-ZSUPLO1)
        ZINF    = ZINFLO1+PDLAT(JROF,JLEV)*(ZINFLO2-ZINFLO1)
        PXF(JROF,JLEV) = ZSUP+PDVER(JROF,JLEV)*(ZINF-ZSUP)
      ENDIF
    ENDDO
  ELSEIF (IQM==2) THEN
    DO JROF=KST,KPROF
      ! compute local min/max
      ZMIN=MIN(ZF111(JROF),ZF211(JROF),ZF121(JROF),ZF221(JROF),&
        &       ZF112(JROF),ZF212(JROF),ZF122(JROF),ZF222(JROF))
      ZMAX=MAX(ZF111(JROF),ZF211(JROF),ZF121(JROF),ZF221(JROF),&
        &       ZF112(JROF),ZF212(JROF),ZF122(JROF),ZF222(JROF))
      ! limit
      PXF(JROF,JLEV)=MAX(ZMIN,MIN(ZMAX,PXF(JROF,JLEV)))
    ENDDO
  ELSEIF ((IQM==3)) THEN
    DO JROF=KST,KPROF
      ZMIN=MIN(ZF111(JROF),ZF211(JROF),ZF121(JROF),ZF221(JROF),&
        &       ZF112(JROF),ZF212(JROF),ZF122(JROF),ZF222(JROF))
      ZMAX=MAX(ZF111(JROF),ZF211(JROF),ZF121(JROF),ZF221(JROF),&
        &       ZF112(JROF),ZF212(JROF),ZF122(JROF),ZF222(JROF))
      ZXF=PXF(JROF,JLEV)+ZSURPL(JROF)
      PXF(JROF,JLEV)=MAX(ZMIN,MIN(ZMAX,ZXF))
      ZDMAX=MAX(ZXF,ZMAX) - ZMAX
      ZDMIN=MIN(ZXF,ZMIN) - ZMIN
      IF(ABS(ZDMAX) >= ABS(ZDMIN)) THEN 
        ZSURPL(JROF)=ZDMAX*ZRVETA(JLEV)
      ELSE
        ZSURPL(JROF)=ZDMIN*ZRVETA(JLEV)
      ENDIF
    ENDDO
  ENDIF
ENDDO

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAITRIQM3D',1,ZHOOK_HANDLE)
END SUBROUTINE LAITRIQM3D
