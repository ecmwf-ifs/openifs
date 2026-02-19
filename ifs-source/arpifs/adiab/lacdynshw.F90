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

SUBROUTINE LACDYNSHW(YDGEM,YDDYN,YDDYNA,CDCONF,KPROMA,KSTART,KPROF,PDT,&
 & PUT9,PVT9,PSPT9,PDIVT9,PSPT9L,PSPT9M,&
 & PUT0,PVT0,PSPT0,PDIVT0,PSPT0L,PSPT0M,&
 & YDGSGEOM,POROG,POROGL,POROGM,&
 & PSPNLT9,&
 & PUL9,PVL9,PSPL9,PUL0,PVL0,PSPL0,PUT1,PVT1,PSPT1,&
 & PURL0,PVRL0,PUSI,PVSI,PURL,PVRL)  

USE YOMGEM   , ONLY : TGEM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : N2DINI
USE YOMCT3   , ONLY : NSTEP
USE YOMDYNA  , ONLY : TDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMGSGEOM, ONLY : TGSGEOM


!**** *LACDYNSHW* - Semi-Lagrangian scheme. - 2-D shallow-water model.
!                   Computation of the t and t-dt useful quantities
!                   at grid-points.

!     Purpose.
!     --------
!          This subroutine computes the equation quantities at t-dt and t
!          at each grid-point of the colocation grid (Gauss grid).

!**   Interface.
!     ----------
!        *CALL* *LACDYNSHW(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          CDCONF  - configuration.
!          KPROMA  - horizontal dimension.
!          KSTART  - first element of work.
!          KPROF   - depth of work.
!          PDT     - For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!          PUT9    - U-component of the wind at t-dt.
!          PVT9    - V-component of the wind at t-dt.
!          PSPT9   - equivalent height "Phi" at t-dt (useless if 2TLSL).
!          PDIVT9  - divergence at t-dt (useless if 2TLSL).
!          PSPT9L  - zonal derivative of equivalent height "Phi" at t-dt.
!          PSPT9M  - meridian derivative of equivalent height "Phi" at t-dt.
!          PUT0    - U-component of the wind at t.
!          PVT0    - V-component of the wind at t.
!          PSPT0   - equivalent height "Phi" at t.
!          PDIVT0  - divergence at t.
!          PSPT0L  - zonal derivative of equivalent height "Phi" at t.
!          PSPT0M  - meridian derivative of equivalent height "Phi" at t.
!          YDGSGEOM- grid point geometry.
!          POROG   - grid-point surface orography
!          POROGL  - zonal component of the orography gradient
!          POROGM  - merid component of the orography gradient

!        INPUT/OUTPUT:
!          PSPNLT9 - nonlinear term in continuity equation at t-dt.
!          PUL9    - SL quantity to be interpolated at O for U-wind eqn.
!          PVL9    - SL quantity to be interpolated at O for V-wind eqn.
!          PSPL9   - SL quantity to be interpolated at O for continuity eqn.
!          PUL0    - second SL quantity to be
!                    interpolated (linearly) at O if KWLAG=3,
!                    for the U-wind equation.
!          PVL0    - second SL quantity to be
!                    interpolated (linearly) at O if KWLAG=3,
!                    for the V-wind equation.
!          PSPL0   - second SL quantity to be
!                    interpolated (linearly) at O if NVLAG=3,
!                    for the continuity equation.
!          PUT1    - t+dt term and other final point terms for U-wind eqn.
!          PVT1    - t+dt term and other final point terms for V-wind eqn.
!          PSPT1   - t+dt term and other final point terms for continuity eqn.

!        OUTPUT:
!          PURL0   - U-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be interpolated)
!          PVRL0   - V-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be interpolated)
!          PUSI    - Part of the t semi-implicit term evaluated at the
!                    final point in the U-wind equation.
!          PVSI    - Part of the t semi-implicit term evaluated at the
!                    final point in the V-wind equation.
!          PURL    - U-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be used at the arrival point)
!          PVRL    - V-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be used at the arrival point)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      None
!      Called by CPG2.

!     Reference.
!     ----------
!        Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD after CPG2 and LACDYN.
!      Original : 98-08-05

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 04-02-03 by C.Temperton (bugfix: set T1 fields to zero)
!      Modified 06-02-06 by C.Temperton (bugfixes for 2TL and orography)
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM)        , INTENT(IN)    :: YDGEM
TYPE(TDYN)         ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)        ,INTENT(IN)    :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROMA 
CHARACTER(LEN=1)  , INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM), INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   , INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   , INTENT(IN)    :: PUT9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PVT9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PSPT9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PDIVT9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PSPT9L(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PSPT9M(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PUT0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PVT0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PSPT0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PDIVT0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PSPT0L(KPROMA) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PSPT0M(KPROMA) 
TYPE(TGSGEOM)     , INTENT(IN)    :: YDGSGEOM
REAL(KIND=JPRB)   , INTENT(IN)    :: POROG(KPROMA)
REAL(KIND=JPRB)   , INTENT(IN)    :: POROGL(KPROMA)
REAL(KIND=JPRB)   , INTENT(IN)    :: POROGM(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSPNLT9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PUL9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PVL9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSPL9(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PUL0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PVL0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSPL0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PUT1(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PVT1(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSPT1(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PURL0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PVRL0(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PUSI(KPROMA) 
REAL(KIND=JPRB)   , INTENT(INOUT) :: PVSI(KPROMA) 
REAL(KIND=JPRB)   , INTENT(OUT)   :: PURL(KPROMA) 
REAL(KIND=JPRB)   , INTENT(OUT)   :: PVRL(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZBDT(KPROMA)
REAL(KIND=JPRB) :: ZGAGT9L(KPROMA),ZGAGT9M(KPROMA)
REAL(KIND=JPRB) :: ZGAGT0L(KPROMA),ZGAGT0M(KPROMA)
REAL(KIND=JPRB) :: ZMOY1U(KPROMA)
REAL(KIND=JPRB) :: ZMOY1V(KPROMA)
REAL(KIND=JPRB) :: ZMOYUNL(KPROMA)
REAL(KIND=JPRB) :: ZMOYVNL(KPROMA)
REAL(KIND=JPRB) :: ZMOY1SP(KPROMA)
REAL(KIND=JPRB) :: ZMOYSPNL(KPROMA)
REAL(KIND=JPRB) :: ZSPNLT0(KPROMA)
REAL(KIND=JPRB) :: ZSPSI(KPROMA)
REAL(KIND=JPRB) :: ZJM(KPROMA)
REAL(KIND=JPRB) :: ZJ0(KPROMA)

REAL(KIND=JPRB) :: ZBETA, ZBT, ZCMSLP, ZDTS2, ZESGM, ZESGP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LACDYNSHW',0,ZHOOK_HANDLE)
ASSOCIATE(BETADT=>YDDYN%BETADT, LADVF=>YDDYN%LADVF, LIMPF=>YDDYN%LIMPF, &
 & LSIDG=>YDDYN%LSIDG, NVLAG=>YDDYN%NVLAG, NWLAG=>YDDYN%NWLAG, &
 & RCMSLP0=>YDDYN%RCMSLP0, SIVP=>YDDYN%SIVP, VESL=>YDDYN%VESL, &
 & LELTRA=>YDDYNA%LELTRA,LSETTLS=>YDDYNA%LSETTLS,LSETTLST=>YDDYNA%LSETTLST, &
 & RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

ZDTS2 = 0.5_JPRB*PDT
ZBETA = 0.0_JPRB
IF (CDCONF == 'A'.OR.CDCONF == 'B') THEN
  IF (YDDYNA%LTWOTL) THEN
    ZBETA = 1.0_JPRB
  ELSE
    ZBETA = BETADT
  ENDIF
ELSE
  CALL ABOR1(' LACDYNSHW: Not existing CDCONF ')
ENDIF
!  Williamson et al advection test
IF(N2DINI == 11) ZBETA=0.0_JPRB
ZBT = ZDTS2*ZBETA

IF (LSIDG) THEN
  ZBDT(KSTART:KPROF)=ZBT*SIVP(1)
ELSE
  ZBDT(KSTART:KPROF)=ZBT*SIVP(1)*RSTRET*RSTRET&
   & /(YDGSGEOM%GM(KSTART:KPROF)*YDGSGEOM%GM(KSTART:KPROF))  
ENDIF

ZESGM=1.0_JPRB-VESL
ZESGP=1.0_JPRB+VESL

ZCMSLP=RCMSLP0

PUT1(KSTART:KPROF)=0.0_JPRB
PVT1(KSTART:KPROF)=0.0_JPRB
PSPT1(KSTART:KPROF)=0.0_JPRB

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE LINEAR TERMS FOR SEMI-IMPLICIT SCHEME.
!              ---------------------------------------------------------

!     * Momentum equation.
!       Add semi-implicit Coriolis terms if required (LIMPF=.T.).

ZGAGT0L(KSTART:KPROF)=PSPT0L(KSTART:KPROF)
ZGAGT0M(KSTART:KPROF)=PSPT0M(KSTART:KPROF)
IF (.NOT.YDDYNA%LTWOTL) THEN
  ZGAGT9L(KSTART:KPROF)=PSPT9L(KSTART:KPROF)
  ZGAGT9M(KSTART:KPROF)=PSPT9M(KSTART:KPROF)
ENDIF
IF (LIMPF) THEN
  ZGAGT0L(KSTART:KPROF)=ZGAGT0L(KSTART:KPROF)&
   & -YDGSGEOM%RCORI(KSTART:KPROF)*PVT0(KSTART:KPROF)  
  ZGAGT0M(KSTART:KPROF)=ZGAGT0M(KSTART:KPROF)&
   & +YDGSGEOM%RCORI(KSTART:KPROF)*PUT0(KSTART:KPROF)  
  IF (.NOT.YDDYNA%LTWOTL) THEN
    ZGAGT9L(KSTART:KPROF)=ZGAGT9L(KSTART:KPROF)&
     & -YDGSGEOM%RCORI(KSTART:KPROF)*PVT9(KSTART:KPROF)  
    ZGAGT9M(KSTART:KPROF)=ZGAGT9M(KSTART:KPROF)&
     & +YDGSGEOM%RCORI(KSTART:KPROF)*PUT9(KSTART:KPROF)  
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE WIND COMPONENTS NECESSARY FOR SL TRAJECTORY.
!              ---------------------------------------------------------------

IF (YDDYNA%LTWOTL) THEN

! * 2TL SL scheme.
  IF (LELTRA) THEN
    PURL0(KSTART:KPROF) = PUT0(KSTART:KPROF)&
     & +0.5_JPRB*PDT*(YDGSGEOM%RCORI(KSTART:KPROF)*PVT0(KSTART:KPROF)&
     & -PSPT0L(KSTART:KPROF))  
    PVRL0(KSTART:KPROF) = PVT0(KSTART:KPROF)&
     & -0.5_JPRB*PDT*(YDGSGEOM%RCORI(KSTART:KPROF)*PUT0(KSTART:KPROF)&
     & +PSPT0M(KSTART:KPROF))  
  ELSE
    IF (NSTEP <= 0) THEN
      PURL0(KSTART:KPROF)=PUT0(KSTART:KPROF)
      PVRL0(KSTART:KPROF)=PVT0(KSTART:KPROF)
    ELSE
      IF (LSETTLST) THEN
        PURL0(KSTART:KPROF)=2.0_JPRB*PUT0(KSTART:KPROF)-PUT9(KSTART:KPROF)
        PVRL0(KSTART:KPROF)=2.0_JPRB*PVT0(KSTART:KPROF)-PVT9(KSTART:KPROF)
      ELSE
        PURL0(KSTART:KPROF)=1.5_JPRB*PUT0(KSTART:KPROF)&
         & -0.5_JPRB*PUT9(KSTART:KPROF)  
        PVRL0(KSTART:KPROF)=1.5_JPRB*PVT0(KSTART:KPROF)&
         & -0.5_JPRB*PVT9(KSTART:KPROF)  
      ENDIF
    ENDIF
  ENDIF
  IF (LSETTLST) THEN
    PURL(KSTART:KPROF)=PUT0(KSTART:KPROF)
    PVRL(KSTART:KPROF)=PVT0(KSTART:KPROF)
  ELSE
    PURL(KSTART:KPROF)=PURL0(KSTART:KPROF)
    PVRL(KSTART:KPROF)=PVRL0(KSTART:KPROF)
  ENDIF

ELSE

! * 3TL SL scheme.
  PURL0(KSTART:KPROF)=PUT0(KSTART:KPROF)
  PVRL0(KSTART:KPROF)=PVT0(KSTART:KPROF)
  PURL(KSTART:KPROF)=PUT0(KSTART:KPROF)
  PVRL(KSTART:KPROF)=PVT0(KSTART:KPROF)

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF THE MOMENTUM EQUATION RIGHT-HAND SIDE TERMS.
!              -----------------------------------------------------------

IF (YDDYNA%LTWOTL) THEN

! * Two-time-level scheme:
!   The non linear terms ZMOY1U+PUSI, ZMOY1V+PVSI are assumed to be
!   zero in this case (the orography is treated as part of the linear term)
!   (BETADT=1, LADVF=T or LIMPF=T), so simplifications
!   are already made. 

  PUSI(KSTART:KPROF)=ZBT*ZGAGT0L(KSTART:KPROF)
  PVSI(KSTART:KPROF)=ZBT*ZGAGT0M(KSTART:KPROF)
  IF (NWLAG == 2) THEN
    PUL9(KSTART:KPROF)=PUL9(KSTART:KPROF)+PUT0(KSTART:KPROF)&
     & -ZESGM*PUSI(KSTART:KPROF)  
    PVL9(KSTART:KPROF)=PVL9(KSTART:KPROF)+PVT0(KSTART:KPROF)&
     & -ZESGM*PVSI(KSTART:KPROF)  
    PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
    PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
  ELSEIF(NWLAG == 3) THEN
    PUL9(KSTART:KPROF)=PUL9(KSTART:KPROF)+PUT0(KSTART:KPROF)
    PVL9(KSTART:KPROF)=PVL9(KSTART:KPROF)+PVT0(KSTART:KPROF)
    PUL0(KSTART:KPROF)=PUL0(KSTART:KPROF)-ZESGM*PUSI(KSTART:KPROF)
    PVL0(KSTART:KPROF)=PVL0(KSTART:KPROF)-ZESGM*PVSI(KSTART:KPROF)
    PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
    PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
  ENDIF

ELSE

! * Three-time-level scheme:
!   ZMOYUNL and ZMOYVNL can be non zero. Note that when BETADT=1,
!   and (LADVF=T or LIMPF=T) these two quantities are zero.

  IF (LADVF) THEN
    ZMOY1U(KSTART:KPROF)=-ZDTS2*ZGAGT0L(KSTART:KPROF)
    ZMOY1V(KSTART:KPROF)=-ZDTS2*ZGAGT0M(KSTART:KPROF)
  ELSE
    ZMOY1U(KSTART:KPROF)=-ZDTS2*(ZGAGT0L(KSTART:KPROF)&
     & -YDGSGEOM%RCORI(KSTART:KPROF)*PVT0(KSTART:KPROF))  
    ZMOY1V(KSTART:KPROF)=-ZDTS2*(ZGAGT0M(KSTART:KPROF)&
     & +YDGSGEOM%RCORI(KSTART:KPROF)*PUT0(KSTART:KPROF))  
  ENDIF
  PUSI(KSTART:KPROF)=ZBT*ZGAGT0L(KSTART:KPROF)
  PVSI(KSTART:KPROF)=ZBT*ZGAGT0M(KSTART:KPROF)
  ZMOYUNL(KSTART:KPROF)=ZMOY1U(KSTART:KPROF)+PUSI(KSTART:KPROF)
  ZMOYVNL(KSTART:KPROF)=ZMOY1V(KSTART:KPROF)+PVSI(KSTART:KPROF)
  IF (NWLAG == 2) THEN
    PUT1(KSTART:KPROF)=PUT1(KSTART:KPROF)+ZESGP*ZMOYUNL(KSTART:KPROF)
    PVT1(KSTART:KPROF)=PVT1(KSTART:KPROF)+ZESGP*ZMOYVNL(KSTART:KPROF)
    PUL9(KSTART:KPROF)=PUL9(KSTART:KPROF)+PUT9(KSTART:KPROF)&
     & +ZESGM*(ZMOYUNL(KSTART:KPROF)-ZBT*ZGAGT9L(KSTART:KPROF))
    PVL9(KSTART:KPROF)=PVL9(KSTART:KPROF)+PVT9(KSTART:KPROF)&
     & +ZESGM*ZMOYVNL(KSTART:KPROF)&
     & +ZESGM*(-ZBT*ZGAGT9M(KSTART:KPROF))  
  ELSEIF(NWLAG == 3) THEN
    PUT1(KSTART:KPROF)=PUT1(KSTART:KPROF)+ZESGP*ZMOYUNL(KSTART:KPROF)
    PVT1(KSTART:KPROF)=PVT1(KSTART:KPROF)+ZESGP*ZMOYVNL(KSTART:KPROF)
    PUL9(KSTART:KPROF)=PUL9(KSTART:KPROF)+PUT9(KSTART:KPROF)
    PVL9(KSTART:KPROF)=PVL9(KSTART:KPROF)+PVT9(KSTART:KPROF)
    PUL0(KSTART:KPROF)=PUL0(KSTART:KPROF)&
     & +ZESGM*ZMOYUNL(KSTART:KPROF)&
     & +ZESGM*(-ZBT*ZGAGT9L(KSTART:KPROF))  
    PVL0(KSTART:KPROF)=PVL0(KSTART:KPROF)&
     & +ZESGM*ZMOYVNL(KSTART:KPROF)&
     & +ZESGM*(-ZBT*ZGAGT9M(KSTART:KPROF))  
  ENDIF
  PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
  PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
ENDIF

IF(LADVF) THEN
  PUT1(KSTART:KPROF)=PUT1(KSTART:KPROF)-YDGSGEOM%GOMVRL(KSTART:KPROF)
  PVT1(KSTART:KPROF)=PVT1(KSTART:KPROF)-YDGSGEOM%GOMVRM(KSTART:KPROF)
ENDIF

!     ------------------------------------------------------------------

!*       5.    COMPUTATION OF THE CONTINUITY EQUATION RIGHT-HAND SIDE TERMS.
!              -------------------------------------------------------------

IF (YDDYNA%LTWOTL) THEN

! * Two-time-level scheme:
!   Is coded only for conventional formulation of continuity eqn (NVLAG>0)

  ZMOY1SP(KSTART:KPROF)=-ZDTS2*PDIVT0(KSTART:KPROF)&
   & *(PSPT0(KSTART:KPROF)-POROG(KSTART:KPROF))&
   & +ZDTS2*ZCMSLP*(PUT0(KSTART:KPROF)*POROGL(KSTART:KPROF)&
   & +PVT0(KSTART:KPROF)*POROGM(KSTART:KPROF))  
  ZSPSI(KSTART:KPROF)=ZBDT(KSTART:KPROF)*PDIVT0(KSTART:KPROF)
  ZSPNLT0(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
  IF (NSTEP <= 0) THEN
    ZMOYSPNL(KSTART:KPROF)=ZSPNLT0(KSTART:KPROF)
  ELSE
    IF (LSETTLS.OR.LELTRA) THEN
      ZMOYSPNL(KSTART:KPROF)=&
       & (1.0_JPRB+ZESGM)*ZSPNLT0(KSTART:KPROF)-PSPNLT9(KSTART:KPROF)
    ELSE
      ZMOYSPNL(KSTART:KPROF)=&
       & 1.5_JPRB*ZSPNLT0(KSTART:KPROF)-0.5_JPRB*PSPNLT9(KSTART:KPROF)
    ENDIF
  ENDIF
  IF (NVLAG == 2) THEN
    PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)&
     & +ZESGP*ZMOYSPNL(KSTART:KPROF)&
     & +(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF)  
    PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)&
     & +(PSPT0(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))&
     & -ZESGM*ZSPSI(KSTART:KPROF)+ZESGM*ZMOYSPNL(KSTART:KPROF)  
  ELSEIF(NVLAG == 3) THEN
    IF (LSETTLS) THEN
      PSPL0(KSTART:KPROF)=PSPL0(KSTART:KPROF)&
       & +ZMOYSPNL(KSTART:KPROF)-ZESGM*ZSPSI(KSTART:KPROF)  
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)&
       & +ZESGP*ZSPNLT0(KSTART:KPROF)&
       & +(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF)
    ELSE
      PSPL0(KSTART:KPROF)=PSPL0(KSTART:KPROF)&
       & +ZESGM*ZMOYSPNL(KSTART:KPROF)-ZESGM*ZSPSI(KSTART:KPROF)  
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)&
       & +ZESGP*ZMOYSPNL(KSTART:KPROF)&
       & +(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF)
    ENDIF
    PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)&
     & +(PSPT0(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))
  ENDIF
  PSPNLT9(KSTART:KPROF)=ZSPNLT0(KSTART:KPROF)

ELSE

! * Three-time-level scheme: 

  IF (NVLAG > 0) THEN

!   * Conventional formulation of continuity equation:

    ZMOY1SP(KSTART:KPROF)=-ZDTS2*PDIVT0(KSTART:KPROF)&
     & *PSPT0(KSTART:KPROF)&
     & +ZDTS2*ZCMSLP*(PUT0(KSTART:KPROF)*POROGL(KSTART:KPROF)&
     & +PVT0(KSTART:KPROF)*POROGM(KSTART:KPROF))  
    ZSPSI(KSTART:KPROF)=ZBDT(KSTART:KPROF)*PDIVT0(KSTART:KPROF)
    ZMOYSPNL(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
    IF (NVLAG == 2) THEN
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)&
       & +ZESGP*ZMOYSPNL(KSTART:KPROF)&
       & +(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF)  
      PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)&
       & +(PSPT9(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))&
       & +ZESGM*ZMOYSPNL(KSTART:KPROF)&
       & +ZESGM*(-ZBDT(KSTART:KPROF)*PDIVT9(KSTART:KPROF))  
    ELSEIF (NVLAG == 3) THEN
      PSPL0(KSTART:KPROF)=PSPL0(KSTART:KPROF)&
       & +ZESGM*ZMOYSPNL(KSTART:KPROF)&
       & +ZESGM*(-ZBDT(KSTART:KPROF)*PDIVT9(KSTART:KPROF))  
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)&
       & +ZESGP*ZMOYSPNL(KSTART:KPROF)&
       & +(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF)  
      PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)&
       & +(PSPT9(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))  
    ENDIF

  ELSE

!   * Lagrangian formulation of continuity equation:

    ZMOY1SP(KSTART:KPROF)=0.0_JPRB
    ZJM(KSTART:KPROF)=(1.0_JPRB-0.5_JPRB*PDT*PDIVT0(KSTART:KPROF))&
     & /(1.0_JPRB+0.5_JPRB*PDT*PDIVT0(KSTART:KPROF))  
    ZJ0(KSTART:KPROF)=1.0_JPRB/(1.0_JPRB+0.5_JPRB*PDT*PDIVT0(KSTART:KPROF))
    ZSPSI(KSTART:KPROF)=ZBDT(KSTART:KPROF)*PDIVT0(KSTART:KPROF)&
     & *ZJ0(KSTART:KPROF)  
    ZMOYSPNL(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
    IF (NVLAG == -2) THEN
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)&
       & +ZESGP*ZMOYSPNL(KSTART:KPROF)+(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF)  
      PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)&
       & +(PSPT9(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))&
       & *ZJM(KSTART:KPROF)+ZESGM*ZMOYSPNL(KSTART:KPROF)&
       & +ZESGM*(-ZBDT(KSTART:KPROF)*PDIVT9(KSTART:KPROF)&
       & *ZJM(KSTART:KPROF))  
    ENDIF

  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LACDYNSHW',1,ZHOOK_HANDLE)
END SUBROUTINE LACDYNSHW
