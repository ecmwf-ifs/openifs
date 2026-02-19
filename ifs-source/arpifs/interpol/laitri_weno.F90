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
SUBROUTINE LAITRI_WENO(YDDYN,KSLB1,KPROMA,KST,KPROF,KFLEV, &
 & KFLDN,KFLDX,KQM, &
 & PDLAT,PCLA,PDLO,PCLO,KL0,KNOWENO,PCW,PVINTW, &
 & PXSL,PXF,PALPHA,PRDETAR,LDQM3DCONS)

! Purpose :
! -------
!   LAITRI_WENO - semi-Lagrangian scheme: tri-dimensional 56-point
!   interpolations (using three overlapping 32-points stencils)
!   with optional quasi-monotonic treatment. Type of high
!   order interpolator fully controlled by weights, low order
!   interpolator always linear. The vertical interpolation 
!   is treated by the WENO technique.

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
!     KNOWENO - special boundary treatment for WENO
!     PCW     - C_k functions for WENO
!     PVINTW  - weights for cubic vertical interpolation
!     PXSL    - quantity to be interpolated
!     PALPHA (OPTIONAL) - ALPHA (known as p) exponent for WENO
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
!   Filip Vana   27-July-2017 after work of A. Craciun, P. Smolikova & J. Masek
!                and the original LAITRI code.
!   (c) ECMWF

! Modifications :
! -------------
!    F. Vana  28-11-2018  passing ALPHA as an argument 
!    F. Vana  24-Jun-2019  QM/QMH fixers
!    F. Vana  29-Jan-2020  (single) precision fix
!    F. Vana   1-Jun-2020  QM3D limiter
!    F. Vana  25-Jan-2021  QM/QMH fixer limited to inner cube

!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : LOPT_RS6K
USE YOMDYN   , ONLY : TDYN

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
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
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCW(KPROMA,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVINTW(KPROMA,KFLEV,9)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KSLB1*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PALPHA
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PRDETAR(KFLEV)
LOGICAL        ,OPTIONAL,INTENT(IN) :: LDQM3DCONS

!------------------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV9L1, IV9L2, IV0L0, IV0L1, IV0L2, IV0L3, IV1L0, IV1L1, IV1L2, IV1L3, &
 & IV2L0, IV2L1, IV2L2, IV2L3, IV3L0, IV3L1, IV3L2, IV3L3, IV4L1, IV4L2
INTEGER(KIND=JPIM) :: ISHIFT1, ISHIFT2
INTEGER(KIND=JPIM) :: JLEV, JROF
INTEGER(KIND=JPIM) :: JJ, J0

REAL(KIND=JPRB) :: ZZF111, ZZF112, ZZF121, ZZF122
REAL(KIND=JPRB) :: ZZF211, ZZF212, ZZF221, ZZF222
REAL(KIND=JPRB) :: ZZF110, ZZF113, ZZF120, ZZF123
REAL(KIND=JPRB) :: ZZF210, ZZF213, ZZF220, ZZF223
REAL(KIND=JPRB) :: ZZ19, ZZ29, ZZ00, ZZ10, ZZ20, ZZ30, ZZ01, ZZ11, ZZ21, ZZ31
REAL(KIND=JPRB) :: ZZ02, ZZ12, ZZ22, ZZ32, ZZ03, ZZ13, ZZ23, ZZ33, ZZ14, ZZ24

REAL(KIND=JPRB) :: ZMIN, ZMAX
REAL(KIND=JPRB) :: Z19(KPROMA), Z29(KPROMA), Z00(KPROMA), Z10(KPROMA), Z20(KPROMA)
REAL(KIND=JPRB) :: Z30(KPROMA), Z01(KPROMA), Z11(KPROMA), Z21(KPROMA), Z31(KPROMA)
REAL(KIND=JPRB) :: Z02(KPROMA), Z12(KPROMA), Z22(KPROMA), Z32(KPROMA), Z03(KPROMA)
REAL(KIND=JPRB) :: Z13(KPROMA), Z23(KPROMA), Z33(KPROMA), Z14(KPROMA), Z24(KPROMA)
REAL(KIND=JPRB) :: Z(KPROMA,0:5)

! Bermejo & Staniforth fixer
REAL(KIND=JPRB) :: ZF(8), ZXF
REAL(KIND=JPRB) :: ZSURPL(KPROMA)
REAL(KIND=JPRB) :: ZRVETA(KFLEV)

! WENO specific variables
REAL(KIND=JPRB)    :: ZEPS, ZSW, ZEPS_MAX
REAL(KIND=JPRB)    :: ZFX(1:3), ZBETA(1:3), ZWW(1:3)
REAL(KIND=JPRB)    :: ZC(1:3,1:10)
REAL(KIND=JPRB)    :: ZALPHA(KFLEV)

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

IF (LHOOK) CALL DR_HOOK('LAITRI_WENO',0,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

! 0. Set up
!    ------

ZEPS=1.E-6_JPRB  ! Note: The value here is not (directly) related to any machine precision
ZEPS_MAX=100._JPRB*EPSILON(ZEPS)

IF (PRESENT(PALPHA)) THEN
  ZALPHA(1:KFLEV)=PALPHA
ELSE
  ZALPHA(1:KFLEV)=YDDYN%RALPHA
ENDIF
! Treatment near the top (to avoid instabilitieas we don't allow any overshoots there)
DO JLEV=1,MIN(KFLEV,YDDYN%NLEV_ZALPHA)  ! done for uppermost NLEV_ZALPHA levels
  ZALPHA(JLEV)=MAX(ZALPHA(JLEV),YDDYN%RALPHA_TOP)
ENDDO

IF (KQM == -1) THEN
  ! Consistency check for Bermejo & Staniforth fixer
  IF (PRESENT(PRDETAR)) THEN
    ZRVETA(1)=1.0_JPRB
    ZRVETA(2:KFLEV)=PRDETAR(2:KFLEV)
  ELSE
    CALL ABOR1(' LAITRI_WENO : MISSING ARGUMENT PRDETAR')
  ENDIF
  IF (.NOT.PRESENT(LDQM3DCONS)) THEN
    CALL ABOR1(' LAITRI_WENO : MISSING ARGUMENT LDQM3DCONS')
  ENDIF
  ! Initialization
  ZSURPL(KST:KPROF)=0.0_JPRB
ENDIF

! beta_k parameters
! Could go to setup, provided the compiler ensures it is used from cache.
SELECT CASE (YDDYN%NEDER)
CASE (1)
  !------------------------------D1+D2+D3-------------------------------
  ZC(1,1:10)=(/ 61.0_JPRB/45.0_JPRB,331.0_JPRB/30.0_JPRB,331.0_JPRB/30.0_JPRB, &
   & 61.0_JPRB/45.0_JPRB,-141.0_JPRB/20.0_JPRB,179.0_JPRB/30.0_JPRB,           &
   & -293.0_JPRB/180.0_JPRB,-1259.0_JPRB/60.0_JPRB,179.0_JPRB/30.0_JPRB,       &
   & -141.0_JPRB/20.0_JPRB /)
  ZC(2,1:10)=(/ 61.0_JPRB/45.0_JPRB,248.0_JPRB/15.0_JPRB,721.0_JPRB/30.0_JPRB, &
   & 407.0_JPRB/90.0_JPRB,-553.0_JPRB/60.0_JPRB,103.0_JPRB/10.0_JPRB,          &
   & -683.0_JPRB/180.0_JPRB,-2309.0_JPRB/60.0_JPRB,439.0_JPRB/30.0_JPRB,       &
   & -1193.0_JPRB/60.0_JPRB /) 
  ZC(3,1:10)=(/ 407.0_JPRB/90.0_JPRB,721.0_JPRB/30.0_JPRB,248.0_JPRB/15.0_JPRB, &
   & 61.0_JPRB/45.0_JPRB,-1193.0_JPRB/60.0_JPRB,439.0_JPRB/30.0_JPRB,           &
   & -683.0_JPRB/180.0_JPRB,-2309.0_JPRB/60.0_JPRB,103.0_JPRB/10.0_JPRB,        &
   & -553.0_JPRB/60.0_JPRB /) 

CASE (2)
  !---------------------------------D2+D3-------------------------------
  ZC(1,1:10)=(/ 4.0_JPRB/3.0_JPRB,10.0_JPRB,10.0_JPRB,4.0_JPRB/3.0_JPRB, &
   &-7.0_JPRB,6.0_JPRB,-5.0_JPRB/3.0_JPRB,-19.0_JPRB,6.0_JPRB,-7.0_JPRB /)
  ZC(2,1:10)=(/ 4.0_JPRB/3.0_JPRB,16.0_JPRB,22.0_JPRB,10.0_JPRB/3.0_JPRB, &
   & -9.0_JPRB,10.0_JPRB,-11.0_JPRB/3.0_JPRB,-37.0_JPRB,14.0_JPRB,-17.0_JPRB /)
  ZC(3,1:10)=(/ 10.0_JPRB/3.0_JPRB,22.0_JPRB,16.0_JPRB,4.0_JPRB/3.0_JPRB, &
   & -17.0_JPRB,14.0_JPRB,-11.0_JPRB/3.0_JPRB,-37.0_JPRB,10.0_JPRB,-9.0_JPRB /)

CASE (3)
  !---------------------------------D2----------------------------------
  ZC(1,1:10)=(/ 1.0_JPRB/3.0_JPRB,1.0_JPRB,1.0_JPRB,1.0_JPRB/3.0_JPRB,-1.0_JPRB,   & 
   & 0.0_JPRB,1.0_JPRB/3.0_JPRB,-1.0_JPRB,0.0_JPRB,-1.0_JPRB /)
  ZC(2,1:10)=(/ 1.0_JPRB/3.0_JPRB,7.0_JPRB,13.0_JPRB,7.0_JPRB/3.0_JPRB,-3.0_JPRB,  &
   & 4.0_JPRB,-5.0_JPRB/3.0_JPRB,-19.0_JPRB,8.0_JPRB,-11.0_JPRB /)
  ZC(3,1:10)=(/ 7.0_JPRB/3.0_JPRB,13.0_JPRB,7.0_JPRB,1.0_JPRB/3.0_JPRB,-11.0_JPRB, & 
   & 8.0_JPRB,-5.0_JPRB/3.0_JPRB,-19.0_JPRB,4.0_JPRB,-3.0_JPRB /)

CASE (4)
  !---------------------------------D3----------------------------------
  ZC(1,1:10)=(/ 1.0_JPRB,9.0_JPRB,9.0_JPRB,1.0_JPRB,-6.0_JPRB, &
   & 6.0_JPRB,-2.0_JPRB,-18.0_JPRB,6.0_JPRB,-6.0_JPRB /)
  ZC(2,1:10)=(/ 1.0_JPRB,9.0_JPRB,9.0_JPRB,1.0_JPRB,-6.0_JPRB, &
   & 6.0_JPRB,-2.0_JPRB,-18.0_JPRB,6.0_JPRB,-6.0_JPRB /)
  ZC(3,1:10)=(/ 1.0_JPRB,9.0_JPRB,9.0_JPRB,1.0_JPRB,-6.0_JPRB, &
   & 6.0_JPRB,-2.0_JPRB,-18.0_JPRB,6.0_JPRB,-6.0_JPRB /)

CASE (5)
  !--------- smoothness indicators computed using undivided differences ---!
  ZC(1,1:10)=(/ 11.0_JPRB/6.0_JPRB,73.0_JPRB/6.0_JPRB,73.0_JPRB/6.0_JPRB, &
   & 11.0_JPRB/6.0_JPRB,-26.0_JPRB/3.0_JPRB,7.0_JPRB,-2.0_JPRB,           &
   & -68.0_JPRB/3.0_JPRB,7.0_JPRB,-26.0_JPRB/3.0_JPRB /)
  ZC(2,1:10)=(/ 11.0_JPRB/6.0_JPRB,73.0_JPRB/6.0_JPRB,73.0_JPRB/6.0_JPRB, &
   & 11.0_JPRB/6.0_JPRB,-26.0_JPRB/3.0_JPRB,7.0_JPRB,-2.0_JPRB,           &
   & -68.0_JPRB/3.0_JPRB,7.0_JPRB,-26.0_JPRB/3.0_JPRB /)
  ZC(3,1:10)=(/ 11.0_JPRB/6.0_JPRB,73.0_JPRB/6.0_JPRB,73.0_JPRB/6.0_JPRB, &
   & 11.0_JPRB/6.0_JPRB,-26.0_JPRB/3.0_JPRB,7.0_JPRB,-2.0_JPRB,           &
   & -68.0_JPRB/3.0_JPRB,7.0_JPRB,-26.0_JPRB/3.0_JPRB /)

CASE DEFAULT
  CALL ABOR1(' LAITRI_WENO: Unrecognized setting for beta_k ')
END SELECT

! offsets for getting values on interpolation stencil
IV9L1=1-KSLB1
IV9L2=2-KSLB1
IV0L0=0
IV0L1=1
IV0L2=2
IV0L3=3
IV1L0=  KSLB1
IV1L1=1+KSLB1
IV1L2=2+KSLB1
IV1L3=3+KSLB1
IV2L0=IV1L0+KSLB1
IV2L1=IV1L1+KSLB1
IV2L2=IV1L2+KSLB1
IV2L3=IV1L3+KSLB1
IV3L0=IV2L0+KSLB1
IV3L1=IV2L1+KSLB1
IV3L2=IV2L2+KSLB1
IV3L3=IV2L3+KSLB1
IV4L1=IV3L1+KSLB1
IV4L2=IV3L2+KSLB1


! 1. Interpolations
!    --------------

#ifdef cray
  IF(0==1)THEN  ! Use Vector code
#else
  IF (LOPT_RS6K) THEN
#endif

! 56-point interpolations
  DO JLEV=1,KFLEV

    ! interpolations in longitude
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ! interpolations in longitude, stencil level 9
      ISHIFT1=MAX(0,KNOWENO(JROF,JLEV))*KSLB1
      Z19(JROF)=PXSL(KL0(JROF,JLEV,1)+IV9L1+ISHIFT1)+PDLO(JROF,JLEV,1) &
       & *(PXSL(KL0(JROF,JLEV,1)+IV9L2+ISHIFT1)-PXSL(KL0(JROF,JLEV,1)+IV9L1+ISHIFT1))
      Z29(JROF)=PXSL(KL0(JROF,JLEV,2)+IV9L1+ISHIFT1)+PDLO(JROF,JLEV,2) &
       & *(PXSL(KL0(JROF,JLEV,2)+IV9L2+ISHIFT1)-PXSL(KL0(JROF,JLEV,2)+IV9L1+ISHIFT1))

      ! interpolations in longitude, stencil level 0
      Z00(JROF)=PXSL(KL0(JROF,JLEV,0)+IV0L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV0L2)-PXSL(KL0(JROF,JLEV,0)+IV0L1))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF110=PXSL(KL0(JROF,JLEV,1)+IV0L1)
      ZZF210=PXSL(KL0(JROF,JLEV,1)+IV0L2)
      Z10(JROF)=PXSL(KL0(JROF,JLEV,1)+IV0L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF110-PXSL(KL0(JROF,JLEV,1)+IV0L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF210-PXSL(KL0(JROF,JLEV,1)+IV0L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV0L3)-PXSL(KL0(JROF,JLEV,1)+IV0L0))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF120=PXSL(KL0(JROF,JLEV,2)+IV0L1)
      ZZF220=PXSL(KL0(JROF,JLEV,2)+IV0L2)
      Z20(JROF)=PXSL(KL0(JROF,JLEV,2)+IV0L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF120-PXSL(KL0(JROF,JLEV,2)+IV0L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF220-PXSL(KL0(JROF,JLEV,2)+IV0L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV0L3)-PXSL(KL0(JROF,JLEV,2)+IV0L0))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      Z30(JROF)=PXSL(KL0(JROF,JLEV,3)+IV0L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV0L2)-PXSL(KL0(JROF,JLEV,3)+IV0L1))
      ! interpolations in longitude, stencil level 1
      Z01(JROF)=PXSL(KL0(JROF,JLEV,0)+IV1L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV1L2)-PXSL(KL0(JROF,JLEV,0)+IV1L1))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF111=PXSL(KL0(JROF,JLEV,1)+IV1L1)
      ZZF211=PXSL(KL0(JROF,JLEV,1)+IV1L2)
      Z11(JROF)=PXSL(KL0(JROF,JLEV,1)+IV1L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF111-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF211-PXSL(KL0(JROF,JLEV,1)+IV1L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV1L3)-PXSL(KL0(JROF,JLEV,1)+IV1L0))
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z11:
        ZMIN=MINJ(ZZF111,ZZF211)
        ZMAX=MAXJ(ZZF111,ZZF211)
        Z11(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z11(JROF)))
      ENDIF
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF121=PXSL(KL0(JROF,JLEV,2)+IV1L1)
      ZZF221=PXSL(KL0(JROF,JLEV,2)+IV1L2)
      Z21(JROF)=PXSL(KL0(JROF,JLEV,2)+IV1L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF121-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF221-PXSL(KL0(JROF,JLEV,2)+IV1L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV1L3)-PXSL(KL0(JROF,JLEV,2)+IV1L0))
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z21:
        ZMIN=MINJ(ZZF121,ZZF221)
        ZMAX=MAXJ(ZZF121,ZZF221)
        Z21(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z21(JROF)))
      ENDIF
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      Z31(JROF)=PXSL(KL0(JROF,JLEV,3)+IV1L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV1L2)-PXSL(KL0(JROF,JLEV,3)+IV1L1))
      ! interpolations in longitude, stencil level 2
      Z02(JROF)=PXSL(KL0(JROF,JLEV,0)+IV2L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV2L2)-PXSL(KL0(JROF,JLEV,0)+IV2L1))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF112=PXSL(KL0(JROF,JLEV,1)+IV2L1)
      ZZF212=PXSL(KL0(JROF,JLEV,1)+IV2L2)
      Z12(JROF)=PXSL(KL0(JROF,JLEV,1)+IV2L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF112-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF212-PXSL(KL0(JROF,JLEV,1)+IV2L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV2L3)-PXSL(KL0(JROF,JLEV,1)+IV2L0))
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z12:
        ZMIN=MINJ(ZZF112,ZZF212)
        ZMAX=MAXJ(ZZF112,ZZF212)
        Z12(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z12(JROF)))
      ENDIF
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF122=PXSL(KL0(JROF,JLEV,2)+IV2L1)
      ZZF222=PXSL(KL0(JROF,JLEV,2)+IV2L2)
      Z22(JROF)=PXSL(KL0(JROF,JLEV,2)+IV2L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF122-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF222-PXSL(KL0(JROF,JLEV,2)+IV2L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV2L3)-PXSL(KL0(JROF,JLEV,2)+IV2L0))
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z22:
        ZMIN=MINJ(ZZF122,ZZF222)
        ZMAX=MAXJ(ZZF122,ZZF222)
        Z22(JROF)=MAXJ(ZMIN,MINJ(ZMAX,Z22(JROF)))
      ENDIF
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      Z32(JROF)=PXSL(KL0(JROF,JLEV,3)+IV2L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV2L2)-PXSL(KL0(JROF,JLEV,3)+IV2L1))
      ! interpolations in longitude, stencil level 3
      Z03(JROF)=PXSL(KL0(JROF,JLEV,0)+IV3L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV3L2)-PXSL(KL0(JROF,JLEV,0)+IV3L1))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF113=PXSL(KL0(JROF,JLEV,1)+IV3L1)
      ZZF213=PXSL(KL0(JROF,JLEV,1)+IV3L2)
      Z13(JROF)=PXSL(KL0(JROF,JLEV,1)+IV3L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF113-PXSL(KL0(JROF,JLEV,1)+IV3L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF213-PXSL(KL0(JROF,JLEV,1)+IV3L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV3L3)-PXSL(KL0(JROF,JLEV,1)+IV3L0))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ZZF123=PXSL(KL0(JROF,JLEV,2)+IV3L1)
      ZZF223=PXSL(KL0(JROF,JLEV,2)+IV3L2)
      Z23(JROF)=PXSL(KL0(JROF,JLEV,2)+IV3L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF123-PXSL(KL0(JROF,JLEV,2)+IV3L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF223-PXSL(KL0(JROF,JLEV,2)+IV3L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV3L3)-PXSL(KL0(JROF,JLEV,2)+IV3L0))
    ENDDO
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      Z33(JROF)=PXSL(KL0(JROF,JLEV,3)+IV3L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV3L2)-PXSL(KL0(JROF,JLEV,3)+IV3L1))
      ! interpolations in longitude, stencil level 4
      ISHIFT2=MIN(0,KNOWENO(JROF,JLEV))*KSLB1
      Z14(JROF)=PXSL(KL0(JROF,JLEV,1)+IV4L1+ISHIFT2)+PDLO(JROF,JLEV,1) &
       & *(PXSL(KL0(JROF,JLEV,1)+IV4L2+ISHIFT2)-PXSL(KL0(JROF,JLEV,1)+IV4L1+ISHIFT2))
      Z24(JROF)=PXSL(KL0(JROF,JLEV,2)+IV4L1+ISHIFT2)+PDLO(JROF,JLEV,2) &
       & *(PXSL(KL0(JROF,JLEV,2)+IV4L2+ISHIFT2)-PXSL(KL0(JROF,JLEV,2)+IV4L1+ISHIFT2))
    ENDDO

    DO JROF=KST,KPROF
      ! interpolations in latitude, stencil levels 9, 0, 1, 2, 3, 4
      Z(JROF,0)=Z19(JROF)+PDLAT(JROF,JLEV)*(Z29(JROF)-Z19(JROF))
      Z(JROF,1)=Z00(JROF)+PCLA(JROF,JLEV,1)*(Z10(JROF)-Z00(JROF)) &
       & +PCLA(JROF,JLEV,2)*(Z20(JROF)-Z00(JROF)) &
       & +PCLA(JROF,JLEV,3)*(Z30(JROF)-Z00(JROF))
      Z(JROF,2)=Z01(JROF)+PCLA(JROF,JLEV,1)*(Z11(JROF)-Z01(JROF)) &
       & +PCLA(JROF,JLEV,2)*(Z21(JROF)-Z01(JROF)) &
       & +PCLA(JROF,JLEV,3)*(Z31(JROF)-Z01(JROF))
      Z(JROF,3)=Z02(JROF)+PCLA(JROF,JLEV,1)*(Z12(JROF)-Z02(JROF)) &
       & +PCLA(JROF,JLEV,2)*(Z22(JROF)-Z02(JROF)) &
       & +PCLA(JROF,JLEV,3)*(Z32(JROF)-Z02(JROF))
      Z(JROF,4)=Z03(JROF)+PCLA(JROF,JLEV,1)*(Z13(JROF)-Z03(JROF)) &
       & +PCLA(JROF,JLEV,2)*(Z23(JROF)-Z03(JROF)) &
       & +PCLA(JROF,JLEV,3)*(Z33(JROF)-Z03(JROF))
      Z(JROF,5)=Z14(JROF)+PDLAT(JROF,JLEV)*(Z24(JROF)-Z14(JROF))

      ! NOTE: For the consistency reason when assembling the final quintic vertical
      !       interpolation the Z(1) and Z(4) MUST NOT BE SECURED by the QM(H) fixer.
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z(2):
        ZMIN=MINJ(Z11(JROF),Z21(JROF))
        ZMAX=MAXJ(Z11(JROF),Z21(JROF))
        Z(JROF,2)=MAXJ(ZMIN,MINJ(ZMAX,Z(JROF,2)))
        ! bound Z(3):
        ZMIN=MINJ(Z12(JROF),Z22(JROF))
        ZMAX=MAXJ(Z12(JROF),Z22(JROF))
        Z(JROF,3)=MAXJ(ZMIN,MINJ(ZMAX,Z(JROF,3)))
      ENDIF
    ENDDO

    ! final interpolation in vertical by WENO scheme
    DO JROF=KST,KPROF

      IF (YDDYN%LWENOBC.OR.(KNOWENO(JROF,JLEV) == 0)) THEN
        DO JJ=1,3
          SELECT CASE (JJ)
            CASE (1)
              J0 = 1
            CASE (2)
              J0 = 2 + MIN(0,KNOWENO(JROF,JLEV))
            CASE (3)
              J0 =     MAX(0,KNOWENO(JROF,JLEV))
          END SELECT

          ZFX(JJ)=Z(JROF,J0) &
           &  + PVINTW(JROF,JLEV,1+3*(JJ-1))*(Z(JROF,J0+1)-Z(JROF,J0)) &
           &  + PVINTW(JROF,JLEV,2+3*(JJ-1))*(Z(JROF,J0+2)-Z(JROF,J0)) &
           &  + PVINTW(JROF,JLEV,3+3*(JJ-1))*(Z(JROF,J0+3)-Z(JROF,J0))

          ZBETA(JJ)=ZC(JJ,1)*Z(JROF,J0+3)*Z(JROF,J0+3) +ZC(JJ,2)*Z(JROF,J0+2)*Z(JROF,J0+2) &
           &      + ZC(JJ,3)*Z(JROF,J0+1)*Z(JROF,J0+1) +ZC(JJ,4)*Z(JROF,J0)*Z(JROF,J0)     &
           &      + ZC(JJ,5)*Z(JROF,J0+3)*Z(JROF,J0+2) +ZC(JJ,6)*Z(JROF,J0+3)*Z(JROF,J0+1) &
           &      + ZC(JJ,7)*Z(JROF,J0+3)*Z(JROF,J0)   +ZC(JJ,8)*Z(JROF,J0+2)*Z(JROF,J0+1) &
           &      + ZC(JJ,9)*Z(JROF,J0+2)*Z(JROF,J0)   +ZC(JJ,10)*Z(JROF,J0+1)*Z(JROF,J0)

          ZWW(JJ)=(ZBETA(JJ)+ZEPS)**ZALPHA(JLEV) ! Denominator
          ZWW(JJ)=SIGN(MAX(ABS(ZWW(JJ)),ZEPS_MAX),ZWW(JJ))
          ZWW(JJ)=PCW(JROF,JLEV,JJ)/ZWW(JJ)
        ENDDO

        ZSW=ZWW(1)+ZWW(2)+ZWW(3)

        ZWW(1)=ZWW(1)/ZSW
        ZWW(2)=ZWW(2)/ZSW
        !ZWW(3)=ZWW(3)/ZSW  not used

        PXF(JROF,JLEV)=ZWW(1)*(ZFX(1)-ZFX(3)) + ZWW(2)*(ZFX(2)-ZFX(3)) + ZFX(3)

      ELSE

        ! no WENO, just standard one interpolation
        PXF(JROF,JLEV)=Z(JROF,1) &
         &  + PVINTW(JROF,JLEV,1)*(Z(JROF,2)-Z(JROF,1)) &
         &  + PVINTW(JROF,JLEV,2)*(Z(JROF,3)-Z(JROF,1)) &
         &  + PVINTW(JROF,JLEV,3)*(Z(JROF,4)-Z(JROF,1))

      ENDIF

      IF (KQM > 1) THEN
        ! bound PXF:
        ZMIN=MINJ(Z(JROF,2),Z(JROF,3))
        ZMAX=MAXJ(Z(JROF,2),Z(JROF,3))
        PXF(JROF,JLEV)=MAXJ(ZMIN,MINJ(ZMAX,PXF(JROF,JLEV)))
      ENDIF
    ENDDO

  ENDDO

ELSE   ! Vector code

  ! Warning: Always make sure that compiler is moving
  !          the IF (KQM) blocks outside the loops.

  ! 56-point interpolations
  DO JLEV=1,KFLEV

    ! interpolations in longitude
!CDIR NODEP
!cdir gthreorder
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      ! interpolations in longitude, stencil level 9
      ISHIFT1=MAX(0,KNOWENO(JROF,JLEV))*KSLB1
      ZZ19=PXSL(KL0(JROF,JLEV,1)+IV9L1+ISHIFT1)+PDLO(JROF,JLEV,1) &
       & *(PXSL(KL0(JROF,JLEV,1)+IV9L2+ISHIFT1)-PXSL(KL0(JROF,JLEV,1)+IV9L1+ISHIFT1))
      ZZ29=PXSL(KL0(JROF,JLEV,2)+IV9L1+ISHIFT1)+PDLO(JROF,JLEV,2) &
       & *(PXSL(KL0(JROF,JLEV,2)+IV9L2+ISHIFT1)-PXSL(KL0(JROF,JLEV,2)+IV9L1+ISHIFT1))
      ! interpolations in longitude, stencil level 0
      ZZ00=PXSL(KL0(JROF,JLEV,0)+IV0L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV0L2)-PXSL(KL0(JROF,JLEV,0)+IV0L1))
      ZZF110=PXSL(KL0(JROF,JLEV,1)+IV0L1)
      ZZF210=PXSL(KL0(JROF,JLEV,1)+IV0L2)
      ZZ10=PXSL(KL0(JROF,JLEV,1)+IV0L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF110-PXSL(KL0(JROF,JLEV,1)+IV0L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF210-PXSL(KL0(JROF,JLEV,1)+IV0L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV0L3)-PXSL(KL0(JROF,JLEV,1)+IV0L0))
      ZZF120=PXSL(KL0(JROF,JLEV,2)+IV0L1)
      ZZF220=PXSL(KL0(JROF,JLEV,2)+IV0L2)
      ZZ20=PXSL(KL0(JROF,JLEV,2)+IV0L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF120-PXSL(KL0(JROF,JLEV,2)+IV0L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF220-PXSL(KL0(JROF,JLEV,2)+IV0L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV0L3)-PXSL(KL0(JROF,JLEV,2)+IV0L0))
      ZZ30=PXSL(KL0(JROF,JLEV,3)+IV0L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV0L2)-PXSL(KL0(JROF,JLEV,3)+IV0L1))
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
      ZZ03=PXSL(KL0(JROF,JLEV,0)+IV3L1)+PDLO(JROF,JLEV,0) &
       & *(PXSL(KL0(JROF,JLEV,0)+IV3L2)-PXSL(KL0(JROF,JLEV,0)+IV3L1))
      ZZF113=PXSL(KL0(JROF,JLEV,1)+IV3L1)
      ZZF213=PXSL(KL0(JROF,JLEV,1)+IV3L2)
      ZZ13=PXSL(KL0(JROF,JLEV,1)+IV3L0) &
       & +PCLO(JROF,JLEV,1,1)*(ZZF113-PXSL(KL0(JROF,JLEV,1)+IV3L0)) &
       & +PCLO(JROF,JLEV,2,1)*(ZZF213-PXSL(KL0(JROF,JLEV,1)+IV3L0)) &
       & +PCLO(JROF,JLEV,3,1)* &
       & (PXSL(KL0(JROF,JLEV,1)+IV3L3)-PXSL(KL0(JROF,JLEV,1)+IV3L0))
      ZZF123=PXSL(KL0(JROF,JLEV,2)+IV3L1)
      ZZF223=PXSL(KL0(JROF,JLEV,2)+IV3L2)
      ZZ23=PXSL(KL0(JROF,JLEV,2)+IV3L0) &
       & +PCLO(JROF,JLEV,1,2)*(ZZF123-PXSL(KL0(JROF,JLEV,2)+IV3L0)) &
       & +PCLO(JROF,JLEV,2,2)*(ZZF223-PXSL(KL0(JROF,JLEV,2)+IV3L0)) &
       & +PCLO(JROF,JLEV,3,2)* &
       & (PXSL(KL0(JROF,JLEV,2)+IV3L3)-PXSL(KL0(JROF,JLEV,2)+IV3L0))
      ZZ33=PXSL(KL0(JROF,JLEV,3)+IV3L1)+PDLO(JROF,JLEV,3) &
       & *(PXSL(KL0(JROF,JLEV,3)+IV3L2)-PXSL(KL0(JROF,JLEV,3)+IV3L1))
      ! interpolations in longitude, stencil level 4
      ISHIFT2=MIN(0,KNOWENO(JROF,JLEV))*KSLB1
      ZZ14=PXSL(KL0(JROF,JLEV,1)+IV4L1+ISHIFT2)+PDLO(JROF,JLEV,1) &
       & *(PXSL(KL0(JROF,JLEV,1)+IV4L2+ISHIFT2)-PXSL(KL0(JROF,JLEV,1)+IV4L1+ISHIFT2))
      ZZ24=PXSL(KL0(JROF,JLEV,2)+IV4L1+ISHIFT2)+PDLO(JROF,JLEV,2) &
       & *(PXSL(KL0(JROF,JLEV,2)+IV4L2+ISHIFT2)-PXSL(KL0(JROF,JLEV,2)+IV4L1+ISHIFT2))
  
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
  
      ! interpolations in latitude, stencil levels 9, 0, 1, 2, 3, 4
      Z(JROF,0)=ZZ19+PDLAT(JROF,JLEV)*(ZZ29-ZZ19)
      Z(JROF,1)=ZZ00+PCLA(JROF,JLEV,1)*(ZZ10-ZZ00) &
       & +PCLA(JROF,JLEV,2)*(ZZ20-ZZ00) &
       & +PCLA(JROF,JLEV,3)*(ZZ30-ZZ00)
      Z(JROF,2)=ZZ01+PCLA(JROF,JLEV,1)*(ZZ11-ZZ01) &
       & +PCLA(JROF,JLEV,2)*(ZZ21-ZZ01) &
       & +PCLA(JROF,JLEV,3)*(ZZ31-ZZ01)
      Z(JROF,3)=ZZ02+PCLA(JROF,JLEV,1)*(ZZ12-ZZ02) &
       & +PCLA(JROF,JLEV,2)*(ZZ22-ZZ02) &
       & +PCLA(JROF,JLEV,3)*(ZZ32-ZZ02)
      Z(JROF,4)=ZZ03+PCLA(JROF,JLEV,1)*(ZZ13-ZZ03) &
       & +PCLA(JROF,JLEV,2)*(ZZ23-ZZ03) &
       & +PCLA(JROF,JLEV,3)*(ZZ33-ZZ03)
      Z(JROF,5)=ZZ14+PDLAT(JROF,JLEV)*(ZZ24-ZZ14)

      ! NOTE: For the consistency reasons when assembling the final quintic vertical
      !       interpolation the Z(1) and Z(4) MUST NOT BE SECURED by the QM(H) fixer.
      IF ((KQM == 1).OR.(KQM == 2)) THEN
        ! bound Z(2):
        ZMIN=MINJ(ZZ11,ZZ21)
        ZMAX=MAXJ(ZZ11,ZZ21)
        Z(JROF,2)=MAXJ(ZMIN,MINJ(ZMAX,Z(JROF,2)))
        ! bound Z(3):
        ZMIN=MINJ(ZZ12,ZZ22)
        ZMAX=MAXJ(ZZ12,ZZ22)
        Z(JROF,3)=MAXJ(ZMIN,MINJ(ZMAX,Z(JROF,3)))
      ENDIF
    ENDDO

      ! final interpolation in vertical by WENO scheme

      ! manually unrolled code to ensure vectorization
      !    (see original JJ loop in the scalar code above)

    DO JROF=KST,KPROF
      !JJ=1,J0=1
      ZFX(1)=Z(JROF,1) &
       &  + PVINTW(JROF,JLEV,1)*(Z(JROF,2)-Z(JROF,1)) &
       &  + PVINTW(JROF,JLEV,2)*(Z(JROF,3)-Z(JROF,1)) &
       &  + PVINTW(JROF,JLEV,3)*(Z(JROF,4)-Z(JROF,1))


      IF (YDDYN%LWENOBC.OR.(KNOWENO(JROF,JLEV) == 0)) THEN
        ! Skip this for the boundary areas 
        ZBETA(1)=ZC(1,1)*Z(JROF,4)*Z(JROF,4) +ZC(1,2)*Z(JROF,3)*Z(JROF,3) &
         &     + ZC(1,3)*Z(JROF,2)*Z(JROF,2) +ZC(1,4)*Z(JROF,1)*Z(JROF,1) &
         &     + ZC(1,5)*Z(JROF,4)*Z(JROF,3) +ZC(1,6)*Z(JROF,4)*Z(JROF,2) &
         &     + ZC(1,7)*Z(JROF,4)*Z(JROF,1) +ZC(1,8)*Z(JROF,3)*Z(JROF,2) &
         &     + ZC(1,9)*Z(JROF,3)*Z(JROF,1) +ZC(1,10)*Z(JROF,2)*Z(JROF,1)
        ZWW(1)=(ZBETA(1)+ZEPS)**ZALPHA(JLEV) ! Denominator
        ZWW(1)=SIGN(MAX(ABS(ZWW(1)),ZEPS_MAX),ZWW(1))
        ZWW(1)=PCW(JROF,JLEV,1)/ZWW(1)
        !JJ=2,J0=2
        J0=2+MIN(0,KNOWENO(JROF,JLEV))
        ZFX(2)=Z(JROF,J0) &
         &  + PVINTW(JROF,JLEV,4)*(Z(JROF,J0+1)-Z(JROF,J0)) &
         &  + PVINTW(JROF,JLEV,5)*(Z(JROF,J0+2)-Z(JROF,J0)) &
         &  + PVINTW(JROF,JLEV,6)*(Z(JROF,J0+3)-Z(JROF,J0))
        ZBETA(2)=ZC(2,1)*Z(JROF,J0+3)*Z(JROF,J0+3) +ZC(2,2)*Z(JROF,J0+2)*Z(JROF,J0+2) &
         &     + ZC(2,3)*Z(JROF,J0+1)*Z(JROF,J0+1) +ZC(2,4)*Z(JROF,J0)*Z(JROF,J0)     &
         &     + ZC(2,5)*Z(JROF,J0+3)*Z(JROF,J0+2) +ZC(2,6)*Z(JROF,J0+3)*Z(JROF,J0+1) &
         &     + ZC(2,7)*Z(JROF,J0+3)*Z(JROF,J0)   +ZC(2,8)*Z(JROF,J0+2)*Z(JROF,J0+1) &
         &     + ZC(2,9)*Z(JROF,J0+2)*Z(JROF,J0)   +ZC(2,10)*Z(JROF,J0+1)*Z(JROF,J0)
        ZWW(2)=(ZBETA(2)+ZEPS)**ZALPHA(JLEV) ! Denominator
        ZWW(2)=SIGN(MAX(ABS(ZWW(2)),ZEPS_MAX),ZWW(2))
        ZWW(2)=PCW(JROF,JLEV,2)/ZWW(2)
        ! JJ=3,J0=0
        J0=MAX(0,KNOWENO(JROF,JLEV))
        ZFX(3)=Z(JROF,J0) &
         &  + PVINTW(JROF,JLEV,7)*(Z(JROF,J0+1)-Z(JROF,J0)) &
         &  + PVINTW(JROF,JLEV,8)*(Z(JROF,J0+2)-Z(JROF,J0)) &
         &  + PVINTW(JROF,JLEV,9)*(Z(JROF,J0+3)-Z(JROF,J0))
        ZBETA(3)=ZC(3,1)*Z(JROF,J0+3)*Z(JROF,J0+3) +ZC(3,2)*Z(JROF,J0+2)*Z(JROF,J0+2) &
         &     + ZC(3,3)*Z(JROF,J0+1)*Z(JROF,J0+1) +ZC(3,4)*Z(JROF,J0)*Z(JROF,J0)     &
         &     + ZC(3,5)*Z(JROF,J0+3)*Z(JROF,J0+2) +ZC(3,6)*Z(JROF,J0+3)*Z(JROF,J0+1) &
         &     + ZC(3,7)*Z(JROF,J0+3)*Z(JROF,J0)   +ZC(3,8)*Z(JROF,J0+2)*Z(JROF,J0+1) &
         &     + ZC(3,9)*Z(JROF,J0+2)*Z(JROF,J0)   +ZC(3,10)*Z(JROF,J0+1)*Z(JROF,J0)
        ZWW(3)=(ZBETA(3)+ZEPS)**ZALPHA(JLEV) ! Denominator
        ZWW(3)=SIGN(MAX(ABS(ZWW(3)),ZEPS_MAX),ZWW(3))
        ZWW(3)=PCW(JROF,JLEV,3)/ZWW(3)
  
        ZSW=ZWW(1)+ZWW(2)+ZWW(3)
  
        ZWW(1)=ZWW(1)/ZSW 
        ZWW(2)=ZWW(2)/ZSW
        !ZWW(3)=ZWW(3)/ZSW  not used

        PXF(JROF,JLEV)=ZWW(1)*(ZFX(1)-ZFX(3)) + ZWW(2)*(ZFX(2)-ZFX(3)) + ZFX(3)

      ELSE
        PXF(JROF,JLEV)=ZFX(1)
      ENDIF

      IF (KQM > 1) THEN
        ! bound PXF:
        ZMIN=MINJ(Z(JROF,2),Z(JROF,3))
        ZMAX=MAXJ(Z(JROF,2),Z(JROF,3))
        PXF(JROF,JLEV)=MAXJ(ZMIN,MINJ(ZMAX,PXF(JROF,JLEV)))
      ENDIF

    ENDDO
  
  ENDDO

ENDIF

!------------------------------------------------------------------------------
! 8 point quasi-monotone correction (Bermejo & Staniforth)
!   with conservation (when activated)

!------------------------------------------------------------------------------

IF (KQM == -1 ) THEN
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

IF (LHOOK) CALL DR_HOOK('LAITRI_WENO',1,ZHOOK_HANDLE)
END SUBROUTINE LAITRI_WENO
