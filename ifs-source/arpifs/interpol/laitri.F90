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
SUBROUTINE LAITRI(KSLB1,KPROMA,KST,KPROF,KFLEV, &
 & KFLDN,KFLDX,KQM, &
 & PDLAT,PCLA,PDLO,PCLO,KL0,PVINTW, &
 & PXSL,PXF,PRDETAR,LDQM3DCONS)

! Purpose :
! -------
!   LAITRI - semi-Lagrangian scheme: tri-dimensional 32-point
!   interpolations (with optional quasi-monotonic treatment). Type of high
!   order interpolator fully controlled by weights, low order
!   interpolator always linear.

! Interface :
! ---------
!   INPUT:
!     KSLB1  - horizontal dimension for grid-point quantities
!     KPROMA  - horizontal dimension for interpolation point quantities
!     KST     - first element of arrays where computations are performed
!     KPROF   - depth of work
!     KFLEV   - vertical dimension
!     KFLDN   - number of the first field
!     KFLDX   - number of the last field
!     KQM     - index of monotonicity
!              -1: Bermejo, Staniforth MWR 1992, Vol. 20.
!               0: not monotonous interpolation
!               1: horizontally quasi-monotonous interpolation
!               2: quasi-monotonous interpolation
!               3: vertically quasi-monotonous interpolation
!     PDLAT   - distance for horizontal linear interpolations in latitude
!     PCLA    - weights for horizontal cubic interpolations in latitude
!     PDLO    - distances for horizontal linear interpolations
!               in longitude (latitude rows 0, 1, 2, 3)
!     PCLO    - weights for horizontal cubic interpolations in longitude 
!               (latitude rows 1, 2)
!     KL0     - indices of the four western points of the 16 point
!               interpolation grid
!     PVINTW  - weights for cubic vertical interpolation
!     PXSL    - quantity to be interpolated
!     PRDETAR (OPTIONAL) - VRDETAR (vertical coordinate related par)
!     LDQM3DCONS (OPTIONAL) - conservativeness for Bermejo & Staniforth QM limiter
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
!   ??-Feb-1992 K. Yessad, after the subroutine LAGINT3 written by
!   Maurice Imbard, Alain Craplet and Michel Rochas,
!   Meteo-France, CNRM/GMAP.

! Modifications :
! -------------
!   01-Oct-2003 M. Hamrud  CY28 Cleaning.
!   30-Jun-2008 J. Masek   High order interpolator fully driven by weights,
!     computation of all weights moved to (E)LASCAW.
!   03-Sep-2008 K. Yessad  Merge QM with not-QM version.
!   R. El Khatib 07-08-2009 Optimisation directive for NEC
!   09-Nov-2009 F. Vana Split to scalar and vector alternatives.
!   26-Jan-2010 R. El Khatib Fix optimisation directives for NEC
!   08-Sep-2014 R. El Khatib Cleaning + use vector code if not rs6k
!   K. Yessad (July 2014): Move some variables.
!   23-Jun-2014 J. Hague: Select "Vector" alternative  for Cray
!   F. Vana   21-Nov-2017: option KQM = 3 (useful for trajectory research in NL)
!   F. Vana    1-Jun-2020: option KQM = -1 (inlined LAQMLIMITER)
! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! arp/ifs dependencies to be solved later.
USE YOMMP0 , ONLY : LOPT_RS6K

!------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KSLB1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KQM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLA(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(KPROMA,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLO(KPROMA,KFLEV,3,2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMA,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVINTW(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KSLB1*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PRDETAR(KFLEV)
LOGICAL        ,OPTIONAL,INTENT(IN) :: LDQM3DCONS

!------------------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV0L1, IV0L2, IV1L0, IV1L1, IV1L2, IV1L3, &
 & IV2L0, IV2L1, IV2L2, IV2L3, IV3L1, IV3L2
INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZZF111, ZZF112, ZZF121, ZZF122
REAL(KIND=JPRB) :: ZZF211, ZZF212, ZZF221, ZZF222
REAL(KIND=JPRB) :: ZZ0, ZZ1, ZZ2, ZZ3
REAL(KIND=JPRB) :: ZZ01, ZZ02, ZZ10, ZZ11, ZZ12, ZZ13
REAL(KIND=JPRB) :: ZZ20, ZZ21, ZZ22, ZZ23, ZZ31, ZZ32

REAL(KIND=JPRB) :: ZF111, ZF112, ZF121, ZF122
REAL(KIND=JPRB) :: ZF211, ZF212, ZF221, ZF222
REAL(KIND=JPRB) :: Z0, Z1, Z2, Z3
REAL(KIND=JPRB) :: ZMIN, ZMAX
REAL(KIND=JPRB) :: Z01(KPROMA), Z02(KPROMA)
REAL(KIND=JPRB) :: Z10(KPROMA), Z11(KPROMA), Z12(KPROMA), Z13(KPROMA)
REAL(KIND=JPRB) :: Z20(KPROMA), Z21(KPROMA), Z22(KPROMA), Z23(KPROMA)
REAL(KIND=JPRB) :: Z31(KPROMA), Z32(KPROMA)

! Bermejo & Staniforth fixer
REAL(KIND=JPRB) :: ZF(8), ZXF
REAL(KIND=JPRB) :: ZSURPL(KPROMA)
REAL(KIND=JPRB) :: ZRVETA(KFLEV)

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

IF (LHOOK) CALL DR_HOOK('LAITRI',0,ZHOOK_HANDLE)

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

#ifdef cray
   IF(0==1)THEN  ! Use Vector code
#else
#ifdef __INTEL_COMPILER
  IF(.TRUE.) THEN
#else
  IF (LOPT_RS6K) THEN
#endif
#endif

! 32-point interpolations
  DO JLEV=1,KFLEV

#ifdef __INTEL_COMPILER
     DO JROF=KST,KPROF
        ! Uniques prefetches for the loops to L2
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,0)+IV1L1),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,0)+IV1L2),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,0)+IV2L1),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,0)+IV2L2),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV0L1),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV0L2),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV1L0),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV1L1),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV1L2),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV2L0),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV2L1),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV2L2),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV2L3),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV3L1),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,1)+IV3L2),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV0L1),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV0L2),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV1L0),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV1L1),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV1L2),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV1L3),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV2L0),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV2L1),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV2L2),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV2L3),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV3L1),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,2)+IV3L2),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,3)+IV1L1),2)
!         CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,3)+IV1L2),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,3)+IV2L1),2)
        CALL MM_PREFETCH(PXSL(KL0(JROF,JLEV,3)+IV2L2),2)
     ENDDO
#endif
     ! First loop
     ! interpolations in longitude
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        ! interpolations in longitude, stencil level 0
        Z10(JROF)=PXSL(KL0(JROF,JLEV,1)+IV0L1)+PDLO(JROF,JLEV,1) &
             & *(PXSL(KL0(JROF,JLEV,1)+IV0L2)-PXSL(KL0(JROF,JLEV,1)+IV0L1))
        Z20(JROF)=PXSL(KL0(JROF,JLEV,2)+IV0L1)+PDLO(JROF,JLEV,2) &
             & *(PXSL(KL0(JROF,JLEV,2)+IV0L2)-PXSL(KL0(JROF,JLEV,2)+IV0L1))
        ! interpolations in longitude, stencil level 1
        Z01(JROF)=PXSL(KL0(JROF,JLEV,0)+IV1L1)+PDLO(JROF,JLEV,0) &
             & *(PXSL(KL0(JROF,JLEV,0)+IV1L2)-PXSL(KL0(JROF,JLEV,0)+IV1L1))
     ENDDO
     ! Second loop
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        ZF111=PXSL(KL0(JROF,JLEV,1)+IV1L1)
        ZF211=PXSL(KL0(JROF,JLEV,1)+IV1L2)
        Z11(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L0) &
             & +PCLO(JROF,JLEV,1,1)*(ZF111-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
             & +PCLO(JROF,JLEV,2,1)*(ZF211-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
             & +PCLO(JROF,JLEV,3,1)* &
             & (PXSL(KL0(JROF,JLEV,1)+IV1L3)-PXSL(KL0(JROF,JLEV,1)+IV1L0))
        IF ((KQM == 1).OR.(KQM == 2)) THEN
           ! bound Z11:
           ZMIN=MINJ(ZF111,ZF211)
           ZMAX=MAXJ(ZF111,ZF211)
           Z11(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z11(JROF)))
        ENDIF
     ENDDO
     ! Third loop
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        ZF121=PXSL(KL0(JROF,JLEV,2)+IV1L1)
        ZF221=PXSL(KL0(JROF,JLEV,2)+IV1L2)
        Z21(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L0) &
             & +PCLO(JROF,JLEV,1,2)*(ZF121-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
             & +PCLO(JROF,JLEV,2,2)*(ZF221-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
             & +PCLO(JROF,JLEV,3,2)* &
             & (PXSL(KL0(JROF,JLEV,2)+IV1L3)-PXSL(KL0(JROF,JLEV,2)+IV1L0))
        IF ((KQM == 1).OR.(KQM == 2)) THEN
           ! bound Z21:
           ZMIN=MINJ(ZF121,ZF221)
           ZMAX=MAXJ(ZF121,ZF221)
           Z21(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z21(JROF)))
        ENDIF
     ENDDO
     ! Fourth loop
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        Z31(JROF)=PXSL(KL0(JROF,JLEV,3)+IV1L1)+PDLO(JROF,JLEV,3) &
             & *(PXSL(KL0(JROF,JLEV,3)+IV1L2)-PXSL(KL0(JROF,JLEV,3)+IV1L1))
        ! interpolations in longitude, stencil level 2
        Z02(JROF)=PXSL(KL0(JROF,JLEV,0)+IV2L1)+PDLO(JROF,JLEV,0) &
             & *(PXSL(KL0(JROF,JLEV,0)+IV2L2)-PXSL(KL0(JROF,JLEV,0)+IV2L1))
     ENDDO
     ! Fifth loop
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        ZF112=PXSL(KL0(JROF,JLEV,1)+IV2L1)
        ZF212=PXSL(KL0(JROF,JLEV,1)+IV2L2)
        Z12(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L0) &
             & +PCLO(JROF,JLEV,1,1)*(ZF112-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
             & +PCLO(JROF,JLEV,2,1)*(ZF212-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
             & +PCLO(JROF,JLEV,3,1)* &
             & (PXSL(KL0(JROF,JLEV,1)+IV2L3)-PXSL(KL0(JROF,JLEV,1)+IV2L0))
        IF ((KQM == 1).OR.(KQM == 2)) THEN
           ! bound Z12:
           ZMIN=MINJ(ZF112,ZF212)
           ZMAX=MAXJ(ZF112,ZF212)
           Z12(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z12(JROF)))
        ENDIF
     ENDDO
     ! Sixth loop
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        ZF122=PXSL(KL0(JROF,JLEV,2)+IV2L1)
        ZF222=PXSL(KL0(JROF,JLEV,2)+IV2L2)
        Z22(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L0) &
             & +PCLO(JROF,JLEV,1,2)*(ZF122-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
             & +PCLO(JROF,JLEV,2,2)*(ZF222-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
             & +PCLO(JROF,JLEV,3,2)* &
             & (PXSL(KL0(JROF,JLEV,2)+IV2L3)-PXSL(KL0(JROF,JLEV,2)+IV2L0))
        IF ((KQM == 1).OR.(KQM == 2)) THEN
           ! bound Z22:
           ZMIN=MINJ(ZF122,ZF222)
           ZMAX=MAXJ(ZF122,ZF222)
           Z22(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z22(JROF)))
        ENDIF
     ENDDO
     ! Seventh loop
     !DIR$ PREFERVECTOR
     DO JROF=KST,KPROF
        Z32(JROF)=PXSL(KL0(JROF,JLEV,3)+IV2L1)+PDLO(JROF,JLEV,3) &
             & *(PXSL(KL0(JROF,JLEV,3)+IV2L2)-PXSL(KL0(JROF,JLEV,3)+IV2L1))
        ! interpolations in longitude, stencil level 3
        Z13(JROF)=PXSL(KL0(JROF,JLEV,1)+IV3L1)+PDLO(JROF,JLEV,1) &
             & *(PXSL(KL0(JROF,JLEV,1)+IV3L2)-PXSL(KL0(JROF,JLEV,1)+IV3L1))
        Z23(JROF)=PXSL(KL0(JROF,JLEV,2)+IV3L1)+PDLO(JROF,JLEV,2) &
             & *(PXSL(KL0(JROF,JLEV,2)+IV3L2)-PXSL(KL0(JROF,JLEV,2)+IV3L1))
     ENDDO
     ! Eigth loop
     DO JROF=KST,KPROF
        ! interpolations in latitude, stencil levels 0, 1, 2, 3
        Z0=Z10(JROF)+PDLAT(JROF,JLEV)*(Z20(JROF)-Z10(JROF))
        Z1=Z01(JROF)+PCLA(JROF,JLEV,1)*(Z11(JROF)-Z01(JROF)) &
             & +PCLA(JROF,JLEV,2)*(Z21(JROF)-Z01(JROF)) &
             & +PCLA(JROF,JLEV,3)*(Z31(JROF)-Z01(JROF))
        Z2=Z02(JROF)+PCLA(JROF,JLEV,1)*(Z12(JROF)-Z02(JROF)) &
             & +PCLA(JROF,JLEV,2)*(Z22(JROF)-Z02(JROF)) &
             & +PCLA(JROF,JLEV,3)*(Z32(JROF)-Z02(JROF))
        Z3=Z13(JROF)+PDLAT(JROF,JLEV)*(Z23(JROF)-Z13(JROF))
        IF ((KQM == 1).OR.(KQM == 2)) THEN
           ! bound Z1:
           ZMIN=MINJ(Z11(JROF),Z21(JROF))
           ZMAX=MAXJ(Z11(JROF),Z21(JROF))
           Z1=MAXJ(ZMIN,MINJ(ZMAX,Z1))
           ! bound Z2:
           ZMIN=MINJ(Z12(JROF),Z22(JROF))
           ZMAX=MAXJ(Z12(JROF),Z22(JROF))
           Z2=MAXJ(ZMIN,MINJ(ZMAX,Z2))
        ENDIF
        ! final interpolation in vertical
        PXF(JROF,JLEV)=Z0 &
             & +PVINTW(JROF,JLEV,1)*(Z1-Z0) &
             & +PVINTW(JROF,JLEV,2)*(Z2-Z0) &
             & +PVINTW(JROF,JLEV,3)*(Z3-Z0)  
        IF (KQM > 1) THEN
           ! bound PXF:
           ZMIN=MINJ(Z1,Z2)
           ZMAX=MAXJ(Z1,Z2)
           PXF(JROF,JLEV)=MAXJ(ZMIN,MINJ(ZMAX,PXF(JROF,JLEV)))
        ENDIF
     ENDDO

  ENDDO

ELSE   ! Vector code

  ! Warning: Always make sure that compiler is moving
  !          the IF (KQM) blocks outside the loops.

  ! 32-point interpolations
  DO JLEV=1,KFLEV

    ! interpolations in longitude
!CDIR NODEP
!cdir gthreorder
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ! interpolations in longitude, stencil level 0
      ZZ10=PXSL(KL0(JROF,JLEV,1)+IV0L1)+PDLO(JROF,JLEV,1) &
       & *(PXSL(KL0(JROF,JLEV,1)+IV0L2)-PXSL(KL0(JROF,JLEV,1)+IV0L1))
      ZZ20=PXSL(KL0(JROF,JLEV,2)+IV0L1)+PDLO(JROF,JLEV,2) &
       & *(PXSL(KL0(JROF,JLEV,2)+IV0L2)-PXSL(KL0(JROF,JLEV,2)+IV0L1))
      ! interpolations in longitude, stencil level 1
      ZZ01=PXSL(KL0(JROF,JLEV,0)+IV1L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV1L2)-PXSL(KL0(JROF,JLEV,0)+IV1L1))
      ZZF111=PXSL(KL0(JROF,JLEV,1)+IV1L1)
      ZZF211=PXSL(KL0(JROF,JLEV,1)+IV1L2)
      ZZ11=PXSL(KL0(JROF,JLEV,1)+IV1L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF111-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF211-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV1L3)-PXSL(KL0(JROF,JLEV,1)+IV1L0))
      ZZF121=PXSL(KL0(JROF,JLEV,2)+IV1L1)
      ZZF221=PXSL(KL0(JROF,JLEV,2)+IV1L2)
      ZZ21=PXSL(KL0(JROF,JLEV,2)+IV1L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF121-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF221-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV1L3)-PXSL(KL0(JROF,JLEV,2)+IV1L0))
      ZZ31=PXSL(KL0(JROF,JLEV,3)+IV1L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV1L2)-PXSL(KL0(JROF,JLEV,3)+IV1L1))
      ! interpolations in longitude, stencil level 2
      ZZ02=PXSL(KL0(JROF,JLEV,0)+IV2L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV2L2)-PXSL(KL0(JROF,JLEV,0)+IV2L1))
      ZZF112=PXSL(KL0(JROF,JLEV,1)+IV2L1)
      ZZF212=PXSL(KL0(JROF,JLEV,1)+IV2L2)
      ZZ12=PXSL(KL0(JROF,JLEV,1)+IV2L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF112-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF212-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV2L3)-PXSL(KL0(JROF,JLEV,1)+IV2L0))
      ZZF122=PXSL(KL0(JROF,JLEV,2)+IV2L1)
      ZZF222=PXSL(KL0(JROF,JLEV,2)+IV2L2)
      ZZ22=PXSL(KL0(JROF,JLEV,2)+IV2L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF122-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF222-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV2L3)-PXSL(KL0(JROF,JLEV,2)+IV2L0))
      ZZ32=PXSL(KL0(JROF,JLEV,3)+IV2L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV2L2)-PXSL(KL0(JROF,JLEV,3)+IV2L1))
      ! interpolations in longitude, stencil level 3
      ZZ13=PXSL(KL0(JROF,JLEV,1)+IV3L1)+PDLO(JROF,JLEV,1) &
       & *(PXSL(KL0(JROF,JLEV,1)+IV3L2)-PXSL(KL0(JROF,JLEV,1)+IV3L1))
      ZZ23=PXSL(KL0(JROF,JLEV,2)+IV3L1)+PDLO(JROF,JLEV,2) &
       & *(PXSL(KL0(JROF,JLEV,2)+IV3L2)-PXSL(KL0(JROF,JLEV,2)+IV3L1))
  
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound ZZ11
        ZMIN=MINJ(ZZF111,ZZF211)
        ZMAX=MAXJ(ZZF111,ZZF211)
        ZZ11=MAXJ(ZMIN,MINJ(ZMAX,ZZ11))
        ! bound ZZ21
        ZMIN=MINJ(ZZF121,ZZF221)
        ZMAX=MAXJ(ZZF121,ZZF221)
        ZZ21=MAXJ(ZMIN,MINJ(ZMAX,ZZ21))
        ! bound ZZ12
        ZMIN=MINJ(ZZF112,ZZF212)
        ZMAX=MAXJ(ZZF112,ZZF212)
        ZZ12=MAXJ(ZMIN,MINJ(ZMAX,ZZ12))
        ! bound ZZ22
        ZMIN=MINJ(ZZF122,ZZF222)
        ZMAX=MAXJ(ZZF122,ZZF222)
        ZZ22=MAXJ(ZMIN,MINJ(ZMAX,ZZ22))
      ENDIF
  
      ! interpolations in latitude, stencil levels 0, 1, 2, 3
      ZZ0=ZZ10+PDLAT(JROF,JLEV)*(ZZ20-ZZ10)
      ZZ1=ZZ01+PCLA(JROF,JLEV,1)*(ZZ11-ZZ01) &
       & +PCLA(JROF,JLEV,2)*(ZZ21-ZZ01) &
       & +PCLA(JROF,JLEV,3)*(ZZ31-ZZ01)
      ZZ2=ZZ02+PCLA(JROF,JLEV,1)*(ZZ12-ZZ02) &
       & +PCLA(JROF,JLEV,2)*(ZZ22-ZZ02) &
       & +PCLA(JROF,JLEV,3)*(ZZ32-ZZ02)
      ZZ3=ZZ13+PDLAT(JROF,JLEV)*(ZZ23-ZZ13)
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z1:
        ZMIN=MINJ(ZZ11,ZZ21)
        ZMAX=MAXJ(ZZ11,ZZ21)
        ZZ1=MAXJ(ZMIN,MINJ(ZMAX,ZZ1))
        ! bound Z2:
        ZMIN=MINJ(ZZ12,ZZ22)
        ZMAX=MAXJ(ZZ12,ZZ22)
        ZZ2=MAXJ(ZMIN,MINJ(ZMAX,ZZ2))
      ENDIF
      ! final interpolation in vertical
      PXF(JROF,JLEV)=ZZ0 &
       & +PVINTW(JROF,JLEV,1)*(ZZ1-ZZ0) &
       & +PVINTW(JROF,JLEV,2)*(ZZ2-ZZ0) &
       & +PVINTW(JROF,JLEV,3)*(ZZ3-ZZ0)  
      IF (KQM > 1) THEN
        ! bound PXF:
        ZMIN=MINJ(ZZ1,ZZ2)
        ZMAX=MAXJ(ZZ1,ZZ2)
        PXF(JROF,JLEV)=MAXJ(ZMIN,MINJ(ZMAX,PXF(JROF,JLEV)))
      ENDIF
    ENDDO
  
  ENDDO

ENDIF
!------------------------------------------------------------------------------

! 8 point quasi-monotone correction (Bermejo & Staniforth) with conservation

IF (KQM == -1) THEN
  ! Consistency check 
  IF (PRESENT(PRDETAR)) THEN
    ZRVETA(1)=1.0_JPRB
    ZRVETA(2:KFLEV)=PRDETAR(2:KFLEV)
  ELSE
    CALL ABOR1(' LAITRI : MISSING ARGUMENT PRDETAR')
  ENDIF
  IF (.NOT.PRESENT(LDQM3DCONS)) THEN
    CALL ABOR1(' LAITRI : MISSING ARGUMENT LDQM3DCONS')
  ENDIF
  ! Initialization
  ZSURPL(KST:KPROF)=0.0_JPRB
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ! get surounding point values
      ZF(1)=PXSL(KL0(JROF,JLEV,1)+IV1L1) ! level above d.p. NW point
      ZF(2)=PXSL(KL0(JROF,JLEV,1)+IV1L2) ! level above d.p. NE point
      ZF(3)=PXSL(KL0(JROF,JLEV,2)+IV1L1) ! level above d.p. SW point
      ZF(4)=PXSL(KL0(JROF,JLEV,2)+IV1L2) ! level above d.p. SE point
      ZF(5)=PXSL(KL0(JROF,JLEV,1)+IV2L1) ! level below d.p. NW point
      ZF(6)=PXSL(KL0(JROF,JLEV,1)+IV2L2) ! level below d.p. NE point
      ZF(7)=PXSL(KL0(JROF,JLEV,2)+IV2L1) ! level below d.p. SW point
      ZF(8)=PXSL(KL0(JROF,JLEV,2)+IV2L2) ! level below d.p. SE point
      ! compute local min/max
      ZMIN=MINVAL(ZF(1:8))
      ZMAX=MAXVAL(ZF(1:8))
      ! limit & store surplus
      ZXF=PXF(JROF,JLEV)
      IF (LDQM3DCONS ) ZXF=ZXF+ZSURPL(JROF)
        PXF(JROF,JLEV)=MAX(ZMIN,MIN(ZMAX,ZXF))
      IF (LDQM3DCONS ) ZSURPL(JROF)=(ZXF-PXF(JROF,JLEV))*ZRVETA(JLEV)
    ENDDO
  ENDDO
ENDIF

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAITRI',1,ZHOOK_HANDLE)
END SUBROUTINE LAITRI
