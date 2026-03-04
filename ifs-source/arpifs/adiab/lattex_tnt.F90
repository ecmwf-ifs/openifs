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

SUBROUTINE LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,KST,KPROF,KXLAG,PESGP,PESGM,PXT9,PMOY1X,PXSI9,&
 & PXSI0,PXT1,PXL0,PXL9,PXLF9,PSIDDHXT1,PSIDDHXL0)

!**** *LATTEX_TNT*   Semi-Lagrangian scheme.
!                Computation of the t and t-dt useful quantities
!                 at grid-points. Equations for tri-dimensional
!                 variables for 3TL scheme

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LATTEX_TNT(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST     - first element of work.
!          KPROF   - depth of work.
!          KXLAG   - type of SL discretisation.
!          PESGP   - (1 + uncentering factor).
!          PESGM   - (1 - uncentering factor).
!          PXT9    - prognostic variable time t-dt
!          PMOY1X  - full nonlinear model at time t [(Delta t/2) "cursive" A]
!          PXSI9   - semi-implicit linear model at time t-dt
!                    [- (Delta t/2) "cursive" B]

!        INPUT/OUTPUT:
!          PXSI0      - semi-implicit linear model at time t
!                       [- (Delta t/2) "cursive" B]
!          PXT1       - t+dt term and other final point terms 
!          PXL0       - second SL quantity to be interpolated (linearly) at O
!                       if KXLAG>3
!          PXL9       - SL quantity to be interpolated at O 
!                       (if NSPLTHOI=1 with diffusive interpolation)
!          PXLF9      - SL quantity to be interpolated at O with high order
!                       interpolation (only meaningfull with NSPLTHOI /= 0)
!          PSIDDHXT1  - semi-implicit  term at final point terms
!          PSIDDHXL0  - second Semi-implicit linear quantity to be interpolated 
!                       (linearly) at O if KXLAG>3

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTEX.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        J. VIVODA (SHMI/CHMI/LACE) - rationalization of LATTEX

!     Modifications.
!     --------------
!     Original : March 2002
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     K. Yessad Aug 2008: rationalisation of dummy arg interfaces + cleanings
!     F. Vana   15-Oct-2009 : optin NSPLTHOI
!     F. Vana   22-Feb-2011 : optin N[X]LAG=4
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN
USE YOMLDDH      , ONLY : TLDDH

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
INTEGER(KIND=JPIM),INTENT(IN)    :: KXLAG
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMOY1X(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSI9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXSI0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXL9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXLF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSIDDHXT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSIDDHXL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)

!     ------------------------------------------------------------------

REAL(KIND=JPRB)    :: ZMOYXNL(YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM) :: JLEV,JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTEX_TNT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NSPLTHOI=>YDDYN%NSPLTHOI, &
 & LRSIDDH=>YDLDDH%LRSIDDH)
!     ------------------------------------------------------------------

DO JLEV=1,NFLEVG

  DO JROF=KST,KPROF
    ZMOYXNL(JROF)=PMOY1X(JROF,JLEV)+PXSI0(JROF,JLEV)
  ENDDO
  IF (KXLAG == 2) THEN
    IF (NSPLTHOI /= 0) THEN
      DO JROF=KST,KPROF
        PXT1(JROF,JLEV)=PXT1(JROF,JLEV)&
         & +PESGP*ZMOYXNL(JROF)-PESGP*PXSI0(JROF,JLEV)
        PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT9(JROF,JLEV)
        PXLF9(JROF,JLEV)=PXLF9(JROF,JLEV)&
         & +PESGM*ZMOYXNL(JROF)-PESGM*PXSI9(JROF,JLEV)
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        PXT1(JROF,JLEV)=PXT1(JROF,JLEV)&
         & +PESGP*ZMOYXNL(JROF)-PESGP*PXSI0(JROF,JLEV)
        PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT9(JROF,JLEV)&
         & +PESGM*ZMOYXNL(JROF)-PESGM*PXSI9(JROF,JLEV)
      ENDDO
    ENDIF
  ELSEIF(KXLAG >= 3) THEN
    DO JROF=KST,KPROF
      PXL0(JROF,JLEV)=PXL0(JROF,JLEV)&
       & +PESGM*ZMOYXNL(JROF)-PESGM*PXSI9(JROF,JLEV)
      PXT1(JROF,JLEV)=PXT1(JROF,JLEV)&
       & +PESGP*ZMOYXNL(JROF)-PESGP*PXSI0(JROF,JLEV)
      PXL9(JROF,JLEV)=PXL9(JROF,JLEV)+PXT9(JROF,JLEV)
    ENDDO
  ENDIF
  PXSI0(KST:KPROF,JLEV)=PESGP*PXSI0(KST:KPROF,JLEV)

ENDDO

!     ------------------------------------------------------------------

IF (LRSIDDH) THEN
  IF(KXLAG >= 3) THEN
    PSIDDHXL0(KST:KPROF,1:NFLEVG)=PSIDDHXL0(KST:KPROF,1:NFLEVG)&
     & +PESGM*PXSI0(KST:KPROF,1:NFLEVG)-PESGM*PXSI9(KST:KPROF,1:NFLEVG)
  ENDIF
  PSIDDHXT1(KST:KPROF,JLEV)=PSIDDHXT1(KST:KPROF,JLEV)+PESGP*PXSI0(KST:KPROF,JLEV)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTEX_TNT',1,ZHOOK_HANDLE)
END SUBROUTINE LATTEX_TNT

