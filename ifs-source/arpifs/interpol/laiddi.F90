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

SUBROUTINE LAIDDI(KPROMA,KPROMB,KST,KPROF,KFLEV, &
 & KFLDN,KFLDX,KQM, &
 & PCLA,PDLO,PCLO,KL0, &
 & PXSL,PXF)

! Purpose :
! -------
!   LAIDDI - semi-Lagrangian scheme: bi-dimensional 12-point
!   interpolations (in horizontal) with optional quasi-monotonic treatment. 
!   Type of high order interpolator fully controlled by weights, 
!   low order interpolator always linear.

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
!     PCLA    - weights for horizontal cubic interpolation in latitude
!     PDLO    - distances for horizontal linear interpolations 
!               in longitude (latitude rows 0, 3)
!     PCLO    - weights for horizontal cubic interpolations in longitude
!               (latitude rows 1, 2)
!     KL0     - indices of the four western points of the 16 point 
!               interpolation grid
!     PXSL    - quantity to be interpolated.
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
!   ??-Feb-1992 K. Yessad, after the subroutine LAGINT2 written by 
!   Maurice Imbard, Alain Craplet and Michel Rochas,
!   Meteo-France, CNRM/GMAP.

! Modifications :
! -------------
!   01-Oct-2003 M. Hamrud  CY28 Cleaning.
!   30-Jun-2008 J. Masek   High order interpolator fully driven by weights,
!     computation of all weights moved to (E)LASCAW.
!   03-Sep-2008 K. Yessad  Merge QM with not-QM version.
! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLA(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(KPROMB,KFLEV,0:3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLO(KPROMB,KFLEV,3,2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,0:3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KPROMA*(KFLDX-KFLDN+1)) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMB,KFLEV) 

!------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZF11(KPROMB), ZF12(KPROMB), ZF21(KPROMB), ZF22(KPROMB)
REAL(KIND=JPRB) :: Z0(KPROMB), Z1(KPROMB), Z2(KPROMB), Z3(KPROMB), ZMIN, ZMAX
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

IF (LHOOK) CALL DR_HOOK('LAIDDI',0,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

! 1. Interpolations
!    --------------

DO JLEV=1,KFLEV

  ! interpolations in longitude, stencil level 0
  DO JROF=KST,KPROF
    Z0(JROF)=PXSL(KL0(JROF,JLEV,0)+1)+PDLO(JROF,JLEV,0) &
     & *(PXSL(KL0(JROF,JLEV,0)+2)-PXSL(KL0(JROF,JLEV,0)+1))
  ENDDO

  ! interpolations in longitude, stencil level 1
  DO JROF=KST,KPROF
    ZF11(JROF)=PXSL(KL0(JROF,JLEV,1)+1)
    ZF21(JROF)=PXSL(KL0(JROF,JLEV,1)+2)
    Z1(JROF)=PXSL(KL0(JROF,JLEV,1)) &
     & +PCLO(JROF,JLEV,1,1)*(ZF11(JROF)-PXSL(KL0(JROF,JLEV,1))) &
     & +PCLO(JROF,JLEV,2,1)*(ZF21(JROF)-PXSL(KL0(JROF,JLEV,1))) &
     & +PCLO(JROF,JLEV,3,1)*(PXSL(KL0(JROF,JLEV,1)+3)-PXSL(KL0(JROF,JLEV,1)))
  ENDDO
  IF (KQM == 1) THEN
    ! bound Z1:
    DO JROF=KST,KPROF
      ZMIN=MINJ(ZF11(JROF),ZF21(JROF))
      ZMAX=MAXJ(ZF11(JROF),ZF21(JROF))
      Z1(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z1(JROF)))
    ENDDO
  ENDIF

  ! interpolations in longitude, stencil level 2
  DO JROF=KST,KPROF
    ZF12(JROF)=PXSL(KL0(JROF,JLEV,2)+1)
    ZF22(JROF)=PXSL(KL0(JROF,JLEV,2)+2)
    Z2(JROF)=PXSL(KL0(JROF,JLEV,2)) &
     & +PCLO(JROF,JLEV,1,2)*(ZF12(JROF)-PXSL(KL0(JROF,JLEV,2))) &
     & +PCLO(JROF,JLEV,2,2)*(ZF22(JROF)-PXSL(KL0(JROF,JLEV,2))) &
     & +PCLO(JROF,JLEV,3,2)*(PXSL(KL0(JROF,JLEV,2)+3)-PXSL(KL0(JROF,JLEV,2)))
  ENDDO
  IF (KQM == 1) THEN
    ! bound Z2:
    DO JROF=KST,KPROF
      ZMIN=MINJ(ZF12(JROF),ZF22(JROF))
      ZMAX=MAXJ(ZF12(JROF),ZF22(JROF))
      Z2(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z2(JROF)))
    ENDDO
  ENDIF

  ! interpolations in longitude, stencil level 3
  DO JROF=KST,KPROF
    Z3(JROF)=PXSL(KL0(JROF,JLEV,3)+1)+PDLO(JROF,JLEV,3) &
     & *(PXSL(KL0(JROF,JLEV,3)+2)-PXSL(KL0(JROF,JLEV,3)+1))
  ENDDO

  ! final interpolation in latitude
  DO JROF=KST,KPROF
    PXF(JROF,JLEV)= Z0(JROF) &
     & +PCLA(JROF,JLEV,1)*(Z1(JROF)-Z0(JROF)) &
     & +PCLA(JROF,JLEV,2)*(Z2(JROF)-Z0(JROF)) &
     & +PCLA(JROF,JLEV,3)*(Z3(JROF)-Z0(JROF))
  ENDDO
  IF (KQM == 1) THEN
    ! bound PXF:
    DO JROF=KST,KPROF
      ZMIN=MINJ(Z1(JROF),Z2(JROF))
      ZMAX=MAXJ(Z1(JROF),Z2(JROF))
      PXF(JROF,JLEV)=MAXJ(ZMIN,MINJ(ZMAX,PXF(JROF,JLEV)))
    ENDDO
  ENDIF

ENDDO

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIDDI',1,ZHOOK_HANDLE)
END SUBROUTINE LAIDDI

