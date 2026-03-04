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

SUBROUTINE LAIHVT(KPROMA,KPROMB,KST,KPROF,KFLEV, &
 & KFLDN,KFLDX,KQM, &
 & PDLAT,PCLA,PDLO,PCLO,KL0,PVDERW,PDZ, &
 & PXSL,PXF)

! Purpose :
! -------
!   LAIHVT - semi-Lagrangian scheme: tri-dimensional 32-point
!   interpolations (with optional quasi-monotonic treatment). Horizontal
!   high order interpolator fully controlled by weights, low order
!   interpolator always linear. Hermite cubic vertical interpolations.

! Interface :
! ---------
!   INPUT:
!     KPROMA  - horizontal dimension for grid-point quantities
!     KPROMB  - horizontal dimension for interpolation point quantities
!     KST     - first element of arrays where computations are performed
!     KPROF   - depth of work
!     KFLEV   - vertical dimension
!     KFLDN   - number of the first field
!     KFLDX   - number of the last field
!     KQM     - index of monotonicity
!               0: not monotonous interpolation
!               1: horizontally quasi-monotonous interpolation
!               2: quasi-monotonous interpolation
!     PDLAT   - distance for horizontal linear interpolations in latitude
!     PCLA    - weights for horizontal cubic interpolations in latitude
!     PDLO    - distances for horizontal linear interpolations
!               in longitude (latitude rows 0, 1, 2, 3)
!     PCLO    - weights for horizontal cubic interpolations in longitude 
!               (latitude rows 1, 2)
!     KL0     - indices of the four western points of the 16 point
!               interpolation grid
!     PVDERW  - weights for computation of vertical derivatives.
!     PDZ     - weights for Hermite cubic vertical interpolation
!     PXSL    - quantity to be interpolated
!   OUTPUT:
!     PXF     - interpolated variable

! Externals :
! ---------
!   None.

! Method :
! ------
!   See documentation.

! Reference :
! ---------

! Author :
! ------
!   K. YESSAD
!   METEO-FRANCE, CNRM/GMAP.

! Modifications :
! -------------
!   Original : AUGUST 1996.
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   F. Vana       16-Sep-2008 Weights driven interpolation.
!   03-Sep-2008 K. Yessad  Merge QM with not-QM version.
! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

!------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KQM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLA(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLO(KPROMB,KFLEV,3,2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDERW(KPROMB,KFLEV,2,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDZ(KPROMB,KFLEV,4)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KPROMA*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMB,KFLEV)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV0L1, IV0L2, IV1L0, IV1L1, IV1L2, IV1L3, &
 & IV2L0, IV2L1, IV2L2, IV2L3, IV3L1, IV3L2
INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZF111(KPROMB), ZF112(KPROMB), ZF121(KPROMB), ZF122(KPROMB)
REAL(KIND=JPRB) :: ZF211(KPROMB), ZF212(KPROMB), ZF221(KPROMB), ZF222(KPROMB)
REAL(KIND=JPRB) :: Z0(KPROMB), Z1(KPROMB), Z2(KPROMB), Z3(KPROMB)
REAL(KIND=JPRB) :: ZMIN, ZMAX
REAL(KIND=JPRB) :: Z01(KPROMB), Z02(KPROMB)
REAL(KIND=JPRB) :: Z10(KPROMB), Z11(KPROMB), Z12(KPROMB), Z13(KPROMB)
REAL(KIND=JPRB) :: Z20(KPROMB), Z21(KPROMB), Z22(KPROMB), Z23(KPROMB)
REAL(KIND=JPRB) :: Z31(KPROMB), Z32(KPROMB)
REAL(KIND=JPRB) :: ZDV1IS, ZDV1SI

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

REAL(KIND=JPRB) :: X, Y, MAXJ, MINJ
#ifdef RS6K
MAXJ(X,Y)=FSEL(X-Y,X,Y)
MINJ(X,Y)=FSEL(X-Y,Y,X)
#else
MAXJ(X,Y)=MAX(X,Y)
MINJ(X,Y)=MIN(X,Y)
#endif

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIHVT',0,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

! 1. Interpolations
!    --------------

! offsets for getting values on interpolation stencil
IV0L1=1
IV0L2=2
IV1L0=  KPROMA
IV1L1=1+KPROMA
IV1L2=2+KPROMA
IV1L3=3+KPROMA
IV2L0=IV1L0+KPROMA
IV2L1=IV1L1+KPROMA
IV2L2=IV1L2+KPROMA
IV2L3=IV1L3+KPROMA
IV3L1=IV2L1+KPROMA
IV3L2=IV2L2+KPROMA

! 32-point interpolations (horizontal interpolations should remain identical
!  to what is done in LAITRI).

!CDIR OUTERUNROLL=4
DO JLEV=1,KFLEV

  ! interpolations in longitude, stencil level 0
!CDIR NODEP
  DO JROF=KST,KPROF
    Z10(JROF)=PXSL(KL0(JROF,JLEV,1)+IV0L1)+PDLO(JROF,JLEV,1) &
     & *(PXSL(KL0(JROF,JLEV,1)+IV0L2)-PXSL(KL0(JROF,JLEV,1)+IV0L1))
    Z20(JROF)=PXSL(KL0(JROF,JLEV,2)+IV0L1)+PDLO(JROF,JLEV,2) &
     & *(PXSL(KL0(JROF,JLEV,2)+IV0L2)-PXSL(KL0(JROF,JLEV,2)+IV0L1))
  ENDDO

  ! interpolations in longitude, stencil level 1
!CDIR NODEP
  DO JROF=KST,KPROF
    Z01(JROF)=PXSL(KL0(JROF,JLEV,0)+IV1L1)+PDLO(JROF,JLEV,0) &
     & *(PXSL(KL0(JROF,JLEV,0)+IV1L2)-PXSL(KL0(JROF,JLEV,0)+IV1L1))
    ZF111(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L1)
    ZF211(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L2)
    Z11(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L0) &
     & +PCLO(JROF,JLEV,1,1)*(ZF111(JROF)-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
     & +PCLO(JROF,JLEV,2,1)*(ZF211(JROF)-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
     & +PCLO(JROF,JLEV,3,1)* &
     & (PXSL(KL0(JROF,JLEV,1)+IV1L3)-PXSL(KL0(JROF,JLEV,1)+IV1L0))
    ZF121(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L1)
    ZF221(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L2)
    Z21(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L0) &
     & +PCLO(JROF,JLEV,1,2)*(ZF121(JROF)-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
     & +PCLO(JROF,JLEV,2,2)*(ZF221(JROF)-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
     & +PCLO(JROF,JLEV,3,2)* &
     & (PXSL(KL0(JROF,JLEV,2)+IV1L3)-PXSL(KL0(JROF,JLEV,2)+IV1L0))
    Z31(JROF)=PXSL(KL0(JROF,JLEV,3)+IV1L1)+PDLO(JROF,JLEV,3) &
     & *(PXSL(KL0(JROF,JLEV,3)+IV1L2)-PXSL(KL0(JROF,JLEV,3)+IV1L1))
  ENDDO
  IF ((KQM == 1).OR.(KQM == 2)) THEN
    DO JROF=KST,KPROF
      ! bound Z11:
      ZMIN=MINJ(ZF111(JROF),ZF211(JROF))
      ZMAX=MAXJ(ZF111(JROF),ZF211(JROF))
      Z11(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z11(JROF)))
      ! bound Z21:
      ZMIN=MINJ(ZF121(JROF),ZF221(JROF))
      ZMAX=MAXJ(ZF121(JROF),ZF221(JROF))
      Z21(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z21(JROF)))
    ENDDO
  ENDIF

  ! interpolations in longitude, stencil level 2
!CDIR NODEP
  DO JROF=KST,KPROF
    Z02(JROF)=PXSL(KL0(JROF,JLEV,0)+IV2L1)+PDLO(JROF,JLEV,0) &
     & *(PXSL(KL0(JROF,JLEV,0)+IV2L2)-PXSL(KL0(JROF,JLEV,0)+IV2L1))
    ZF112(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L1)
    ZF212(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L2)
    Z12(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L0) &
     & +PCLO(JROF,JLEV,1,1)*(ZF112(JROF)-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
     & +PCLO(JROF,JLEV,2,1)*(ZF212(JROF)-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
     & +PCLO(JROF,JLEV,3,1)* &
     & (PXSL(KL0(JROF,JLEV,1)+IV2L3)-PXSL(KL0(JROF,JLEV,1)+IV2L0))
    ZF122(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L1)
    ZF222(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L2)
    Z22(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L0) &
     & +PCLO(JROF,JLEV,1,2)*(ZF122(JROF)-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
     & +PCLO(JROF,JLEV,2,2)*(ZF222(JROF)-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
     & +PCLO(JROF,JLEV,3,2)* &
     & (PXSL(KL0(JROF,JLEV,2)+IV2L3)-PXSL(KL0(JROF,JLEV,2)+IV2L0))
    Z32(JROF)=PXSL(KL0(JROF,JLEV,3)+IV2L1)+PDLO(JROF,JLEV,3) &
     & *(PXSL(KL0(JROF,JLEV,3)+IV2L2)-PXSL(KL0(JROF,JLEV,3)+IV2L1))
  ENDDO
  IF ((KQM == 1).OR.(KQM == 2)) THEN
    DO JROF=KST,KPROF
      ! bound Z12:
      ZMIN=MINJ(ZF112(JROF),ZF212(JROF))
      ZMAX=MAXJ(ZF112(JROF),ZF212(JROF))
      Z12(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z12(JROF)))
      ! bound Z22:
      ZMIN=MINJ(ZF122(JROF),ZF222(JROF))
      ZMAX=MAXJ(ZF122(JROF),ZF222(JROF))
      Z22(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z22(JROF)))
    ENDDO
  ENDIF

  ! interpolations in longitude, stencil level 3
!CDIR NODEP
  DO JROF=KST,KPROF
    Z13(JROF)=PXSL(KL0(JROF,JLEV,1)+IV3L1)+PDLO(JROF,JLEV,1) &
     & *(PXSL(KL0(JROF,JLEV,1)+IV3L2)-PXSL(KL0(JROF,JLEV,1)+IV3L1))
    Z23(JROF)=PXSL(KL0(JROF,JLEV,2)+IV3L1)+PDLO(JROF,JLEV,2) &
     & *(PXSL(KL0(JROF,JLEV,2)+IV3L2)-PXSL(KL0(JROF,JLEV,2)+IV3L1))
  ENDDO

  ! interpolations in latitude, stencil levels 0, 1, 2, 3
  DO JROF=KST,KPROF
    Z0(JROF)=Z10(JROF)+PDLAT(JROF,JLEV)*(Z20(JROF)-Z10(JROF))
    Z1(JROF)=Z01(JROF)+PCLA(JROF,JLEV,1)*(Z11(JROF)-Z01(JROF)) &
     & +PCLA(JROF,JLEV,2)*(Z21(JROF)-Z01(JROF)) &
     & +PCLA(JROF,JLEV,3)*(Z31(JROF)-Z01(JROF))
    Z2(JROF)=Z02(JROF)+PCLA(JROF,JLEV,1)*(Z12(JROF)-Z02(JROF)) &
     & +PCLA(JROF,JLEV,2)*(Z22(JROF)-Z02(JROF)) &
     & +PCLA(JROF,JLEV,3)*(Z32(JROF)-Z02(JROF))
    Z3(JROF)=Z13(JROF)+PDLAT(JROF,JLEV)*(Z23(JROF)-Z13(JROF))
  ENDDO
  IF ((KQM == 1).OR.(KQM == 2)) THEN
    DO JROF=KST,KPROF
      ! bound Z1:
      ZMIN=MINJ(Z11(JROF),Z21(JROF))
      ZMAX=MAXJ(Z11(JROF),Z21(JROF))
      Z1(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z1(JROF)))
      ! bound Z2:
      ZMIN=MINJ(Z12(JROF),Z22(JROF))
      ZMAX=MAXJ(Z12(JROF),Z22(JROF))
      Z2(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z2(JROF)))
    ENDDO
  ENDIF

  ! final interpolation in vertical
  DO JROF=KST,KPROF
    ZDV1SI=PVDERW(JROF,JLEV,1,1)*(Z1(JROF)-Z0(JROF)) &
     & +PVDERW(JROF,JLEV,2,1)*(Z2(JROF)-Z1(JROF))
    ZDV1IS=PVDERW(JROF,JLEV,1,2)*(Z2(JROF)-Z1(JROF)) &
     & +PVDERW(JROF,JLEV,2,2)*(Z3(JROF)-Z2(JROF))
    PXF(JROF,JLEV)= PDZ(JROF,JLEV,1)*Z1(JROF) &
     & +PDZ(JROF,JLEV,2)*Z2(JROF) &
     & +PDZ(JROF,JLEV,3)*ZDV1SI &
     & +PDZ(JROF,JLEV,4)*ZDV1IS
  ENDDO
  IF (KQM == 2) THEN
    ! bound PXF:
    DO JROF=KST,KPROF
      ZMIN=MINJ(Z1(JROF),Z2(JROF))
      ZMAX=MAXJ(Z1(JROF),Z2(JROF))
      PXF(JROF,JLEV)=MAXJ(ZMIN,MINJ(ZMAX,PXF(JROF,JLEV)))
    ENDDO
  ENDIF

ENDDO

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIHVT',1,ZHOOK_HANDLE)
END SUBROUTINE LAIHVT
