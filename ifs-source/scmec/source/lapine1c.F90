! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LAPINE1C(YDGEOMETRY,YDMODEL,KFLEV,PDT,PETAP&
          &,PQT9,PTKET9,PTT9,PUT9,PVT9,PAT9,PLT9,PIT9&
          &,PKAP,PRESF,PVVEL0,PRCORI,PUG0,PVG0&
          &,PTT0 , PUT0 , PVT0,PQT0,PTKET0,PAT0,PLT0,PIT0&
          &,PTADV , PQADV,PTKEADV , PUADV , PVADV&
          &,PUTEND,PVTEND,PTTEND,PQTEND,PTKETEND,PATEND,PLTEND,PITEND&
          &,PSAVTEND,PGFLSLP,PSLPHY9)

#ifdef DOC
!**** *LAPINE1C* - semi-LAgrangian scheme:
!                  Interface subroutine for interpolations and lagrangian
!                  trends. 

!     Purpose.
!     --------
!           Computes the medium point then interpolations at this point.
!           Computes the origin point then interpolations at this point.
!           Computes semi-lagrangian tendencies and t+dt values.

!**   Interface.
!     ----------
!        *CALL* *LAPINE1C

!        Explicit arguments :
!        --------------------
!        INPUT:
!          KFLEV - Number of levels. 
!          PDT   - Time step.
!          PETAP - d(eta)/dt.
!          PUT9  - Quantity to interpolate at the origin point in the
!                  U-component of the wind equation, (U at t-dt).
!          PVT9  - Quantity to interpolate at the origin point in the
!                  V-component of the wind equation, (V at t-dt).
!          PTT9  - Quantity to interpolate at the origin point in the
!                  temperature equation, (T at t-dt).
!          PQT9  - Quantity to interpolate at the origin point in the
!                  humidity equation, (Q at t-dt).
!          PAT9  - A (cloud fraction) at t-dt.
!          PLT9  - L (liq.water) at t-dt.
!          PIT9  - I (ice) at t-dt.
!          PKAP  - K.
!          PRESF - pressure at full-levels.
!          PVVEL0- dp/dt.
!          PRCORI- Coriolis parameter.
!          PUG0  - Geostrophic wind, U component at time t.
!          PVG0  - Geostrophic wind, V component at time t.
!          PTT0  - T at time t.
!          PUT0  - U at time t.
!          PVT0  - V at time t.
!          PQT0  - Q at time t.
!          PTKET0  - TKE at time t.
!          PAT0  - A at time t.
!          PLT0  - QL at time t.
!          PIT0  - QI at time t.
!          PTADV - T horizontal advection, at time t. 
!          PQADV - Q horizontal advection, at time t.
!          PUADV - U horizontal advection, at time t.
!          PVADV - V horizontal advection, at time t.
!          PSAVTEND - tendencies for GMV variables from previous timestep physics
!          PGFLSLP  - tendencies for GFL variables from previous timestep physics
!        OUTPUT:
!          PUTEND- U-wind tendency.
!          PVTEND- V-wind tendency.
!          PTTEND- T tendency.
!          PQTEND- Q tendency.
!          PATEND- A tendency.
!          PLTEND- L tendency.
!          PITEND- I tendency.
!          PSLPHY9 - interpolated phys tendencies
!        LOCAL:  
!          ZUT1  - U at time t+dt. 
!          ZVT1  - V at time t+dt.
!          ZTT1  - T at time t+dt.
!          ZQT1  - Q at time t+dt.
!         TERMS TO BE INTERPOLATED AT THE MEDIUM POINT:
!          ZT2T0 - Tend. from ( hor.adv. of T + subsidence ).
!          ZQ2Q0 - Tend. from hor.adv. of Q.
!          ZU2U0 - Tend. from (hor.adv. + geost. term) for U.
!          ZV2V0 - Tend. from (hor.adv. + geost. term) for V. 
!         TERMS THAT RESULT FROM THE INTERPOLATION AT THE MEDIUM POINT:
!          ZTL0  - From the interpolation of ZT2T0.
!          ZUL0  - From the interpolation of ZU2U0.
!          ZQL0  - From the interpolation of ZQ2Q0.
!          ZVL0  - From the interpolation of ZV2V0.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.    LARMES1C - Computes the medium point coordinates then
!     ----------               interpolations at this point.
!                   LAINOR1C - Computes the origin point coordinates then
!                              interpolations at this point.
!                   Called by CPG1C

!     Reference.
!     ----------


!     Author.
!     -------
!        Joao Teixeira  *ECMWF*        

!     Modifications.
!     --------------
!        Original    18-05-94
!        J.Teixeira    May-1995  introduction of cloud variables
!                                and two time-level scheme.
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!        F. Vana       10-Jul-2014 Rewritten to follow more recent version of
!                                  model SL scheme  + adding SLPHYS
!     ------------------------------------------------------------------
#endif

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMLOG1C , ONLY : LWADVCLD
USE YOMCT3   , ONLY : NSTEP
USE YOMGPD1C , ONLY : RHST9, RHSQ9,RHSTKE9, RHSU9, RHSV9, WRL9,&
                   &  RHSA9,RHSL9,RHSI9
USE YOMPHYDER, ONLY : MODEL_STATE_TYPE

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM) :: KFLEV

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PDT
REAL(KIND=JPRB) :: PRCORI

!     INPUT
REAL(KIND=JPRB) :: PUT9 (YDGEOMETRY%YRDIMV%NFLEVG) , PVT9  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTT9 (YDGEOMETRY%YRDIMV%NFLEVG) , PQT9  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTKET9 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PAT9 (YDGEOMETRY%YRDIMV%NFLEVG) , PLT9  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PIT9 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PUT0 (YDGEOMETRY%YRDIMV%NFLEVG) , PVT0  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTT0 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PQT0 (YDGEOMETRY%YRDIMV%NFLEVG) , PTKET0  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PAT0 (YDGEOMETRY%YRDIMV%NFLEVG) , PLT0  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PIT0 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PUG0 (YDGEOMETRY%YRDIMV%NFLEVG) , PVG0  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PETAP(YDGEOMETRY%YRDIMV%NFLEVG) , PKAP  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PRESF(YDGEOMETRY%YRDIMV%NFLEVG) , PVVEL0(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTADV(YDGEOMETRY%YRDIMV%NFLEVG) , PQADV (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTKEADV(YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: PUADV(YDGEOMETRY%YRDIMV%NFLEVG) , PVADV (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PSAVTEND(YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND)
REAL(KIND=JPRB) :: PGFLSLP(YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMSLP)

!     OUTPUT
REAL(KIND=JPRB) :: PUTEND(YDGEOMETRY%YRDIMV%NFLEVG), PVTEND(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTTEND(YDGEOMETRY%YRDIMV%NFLEVG), PQTEND(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTKETEND(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PATEND(YDGEOMETRY%YRDIMV%NFLEVG), PLTEND(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PITEND(YDGEOMETRY%YRDIMV%NFLEVG)
TYPE (MODEL_STATE_TYPE)   :: PSLPHY9

!     LOCAL
REAL(KIND=JPRB) :: ZT2T0(YDGEOMETRY%YRDIMV%NFLEVG) , ZQ2Q0 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTKE2Q0(YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZA2A0(YDGEOMETRY%YRDIMV%NFLEVG) , ZL2L0 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZI2I0(YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZU2U0(YDGEOMETRY%YRDIMV%NFLEVG) , ZV2V0 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUT1 (YDGEOMETRY%YRDIMV%NFLEVG) , ZVT1  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT1 (YDGEOMETRY%YRDIMV%NFLEVG) , ZQT1  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTKET1 (YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZAT1 (YDGEOMETRY%YRDIMV%NFLEVG) , ZLT1  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZIT1 (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWF  (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTL0 (YDGEOMETRY%YRDIMV%NFLEVG) , ZUL0(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZQL0 (YDGEOMETRY%YRDIMV%NFLEVG) , ZVL0(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTKEL0 (YDGEOMETRY%YRDIMV%NFLEVG)  
REAL(KIND=JPRB) :: ZAL0 (YDGEOMETRY%YRDIMV%NFLEVG) , ZLL0(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZIL0 (YDGEOMETRY%YRDIMV%NFLEVG)  
REAL(KIND=JPRB) :: ZUPT1(YDGEOMETRY%YRDIMV%NFLEVG) , ZVPT1(YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZLEV(YDGEOMETRY%YRDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: J, JLEV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZALFA, ZBETA, ZCONST1, ZCONST2, ZDTCORI, ZBETAT, ZTPT1, ZB, ZA,&
        & ZQPT1,ZTKEPT1,ZAPT1,ZLPT1,ZIPT1

!     ------------------------------------------------------------------
#include "larcin1c.intfb.h"
#include "larmes1c.intfb.h"
!     ------------------------------------------------------------------
ASSOCIATE(YRSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY, YGFL=>YDMODEL%YRML_GCONF%YGFL, &
& LDSLPHY=>YDMODEL%YRML_PHY_EC%YREPHY%LSLPHY)

!*    1. COMPUTATION OF THE TERMS TO BE INTERPOLATED AT THE MEDIUM 
!        POINT AND AT TIME t: 1) TENDENCY FROM THE HORIZONTAL
!        ADVECTION TERM FOR Q; 2) TEND. FROM HOR.ADV. + SUBSIDENCE
!        TERM FOR T; 3) TEND. FROM HOR.ADV. + GEOST. TERM FOR U AND V.
!        -------------------------------------------------------------


IF (YDMODEL%YRML_DYN%YRDYNA%LTWOTL) THEN

  ZBETA=0.5_JPRB   ! implicit factor for momentum equations
  ZBETAT=0.5_JPRB  ! implicit factor for temperature equation
  ZDTCORI=PDT*PRCORI
  ZA   =ZDTCORI*(1._JPRB-ZBETA)
  ZB   =ZDTCORI*ZBETA
  ZALFA=1._JPRB/(1._JPRB + ZB*ZB)

  DO JLEV=1,KFLEV

    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZUPT1(JLEV)=PUT0(JLEV)
      ZVPT1(JLEV)=PVT0(JLEV)
    ELSE
      ! centered scheme when ZBETA=0.5
      ZUPT1(JLEV)=ZALFA*( (1._JPRB-ZA*ZB)*PUT0(JLEV) + ZDTCORI*PVT0(JLEV)&
        & + ZDTCORI*(ZB*PUG0(JLEV)-PVG0(JLEV)) )
      ZVPT1(JLEV)=PVT0(JLEV) - ZA*PUT0(JLEV) - ZB*ZUPT1(JLEV)&
        & + ZDTCORI*PUG0(JLEV)
    ENDIF
    ZU2U0(JLEV)=PDT*(PUADV(JLEV) + (ZVPT1(JLEV)- PVG0(JLEV))*PRCORI)
    ZV2V0(JLEV)=PDT*(PVADV(JLEV) + (PUG0(JLEV) -ZUPT1(JLEV))*PRCORI)

    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZTPT1=PTT0(JLEV)
    ELSE
      ZCONST1=(PDT*(1._JPRB-ZBETAT)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(PDT*ZBETAT          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(1.0_JPRB+ZCONST1)/(1.0_JPRB-ZCONST2)
      ZTPT1=ZCONST2*PTT0(JLEV)
    ENDIF
    ZT2T0(JLEV)=PDT*(PTADV(JLEV)+(ZTPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV))

    !ZQ2Q0(JLEV)=PDT*PQADV(JLEV)
    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZQPT1=PQT0(JLEV)
    ELSE
      ZCONST1=(PDT*(1._JPRB-ZBETAT)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(PDT*ZBETAT          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(1.0_JPRB+ZCONST1)/(1.0_JPRB-ZCONST2)
      ZQPT1=ZCONST2*PQT0(JLEV)
    ENDIF
    ZQ2Q0(JLEV)=PDT*(PQADV(JLEV)+(ZQPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV))

    !ZTKE2Q0(JLEV)=PDT*PTKEADV(JLEV)
    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZTKEPT1=PTKET0(JLEV)
    ELSE
      ZCONST1=(PDT*(1._JPRB-ZBETAT)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(PDT*ZBETAT          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(1.0_JPRB+ZCONST1)/(1.0_JPRB-ZCONST2)
      ZTKEPT1=ZCONST2*PTKET0(JLEV)
    ENDIF
    ZTKE2Q0(JLEV)=PDT*(PTKEADV(JLEV)+(ZTKEPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV))

    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZAPT1=PAT0(JLEV)
    ELSE
      ZCONST1=(PDT*(1._JPRB-ZBETAT)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(PDT*ZBETAT          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(1.0_JPRB+ZCONST1)/(1.0_JPRB-ZCONST2)
      ZAPT1=ZCONST2*PAT0(JLEV)
    ENDIF
    ZA2A0(JLEV)=PDT*((ZAPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV))

    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZLPT1=PLT0(JLEV)
    ELSE
      ZCONST1=(PDT*(1._JPRB-ZBETAT)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(PDT*ZBETAT          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(1.0_JPRB+ZCONST1)/(1.0_JPRB-ZCONST2)
      ZLPT1=ZCONST2*PLT0(JLEV)
    ENDIF
    ZL2L0(JLEV)=PDT*((ZLPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV))

    IF (YDMODEL%YRML_DYN%YRDYNA%LSETTLS) THEN
      ZIPT1=PIT0(JLEV)
    ELSE
      ZCONST1=(PDT*(1._JPRB-ZBETAT)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(PDT*ZBETAT          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
      ZCONST2=(1.0_JPRB+ZCONST1)/(1.0_JPRB-ZCONST2)
      ZIPT1=ZCONST2*PIT0(JLEV)
    ENDIF
    ZI2I0(JLEV)=PDT*((ZIPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV))
    
    IF ((NSTEP > 0).AND.YDMODEL%YRML_DYN%YRDYNA%LSETTLS ) THEN
      ! extrapolation scheme (settls like)
      ZU2U0(JLEV)=ZU2U0(JLEV)-RHSU9(JLEV)
      ZV2V0(JLEV)=ZV2V0(JLEV)-RHSV9(JLEV)
      ZT2T0(JLEV)=ZT2T0(JLEV)-RHST9(JLEV)
      ZQ2Q0(JLEV)=ZQ2Q0(JLEV)-RHSQ9(JLEV)
      ZTKE2Q0(JLEV)=ZTKE2Q0(JLEV)-RHSTKE9(JLEV)
      ZA2A0(JLEV)=ZA2A0(JLEV)-RHSA9(JLEV)
      ZL2L0(JLEV)=ZL2L0(JLEV)-RHSL9(JLEV)
      ZI2I0(JLEV)=ZI2I0(JLEV)-RHSI9(JLEV)
    ELSE
      ! simple averaging along trajectory
      ZU2U0(JLEV)=0.5_JPRB*ZU2U0(JLEV)
      ZV2V0(JLEV)=0.5_JPRB*ZV2V0(JLEV)
      ZT2T0(JLEV)=0.5_JPRB*ZT2T0(JLEV)
      ZQ2Q0(JLEV)=0.5_JPRB*ZQ2Q0(JLEV)
      ZTKE2Q0(JLEV)=0.5_JPRB*ZTKE2Q0(JLEV)
      ZA2A0(JLEV)=0.5_JPRB*ZA2A0(JLEV)
      ZL2L0(JLEV)=0.5_JPRB*ZL2L0(JLEV)
      ZI2I0(JLEV)=0.5_JPRB*ZI2I0(JLEV)
    ENDIF
    ! storing quantities for next timestep and/or averaging along trajectory
    RHSU9(JLEV)=0.5_JPRB*PDT*(PUADV(JLEV) + (ZVPT1(JLEV)- PVG0(JLEV))*PRCORI)
    RHSV9(JLEV)=0.5_JPRB*PDT*(PVADV(JLEV) + (PUG0(JLEV) -ZUPT1(JLEV))*PRCORI)
    RHST9(JLEV)=0.5_JPRB*PDT*(PTADV(JLEV) +&
       & (ZTPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV) )
    RHSQ9(JLEV)=0.5_JPRB*PDT*(PQADV(JLEV) +&
       & (ZQPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV) )
    RHSTKE9(JLEV)=0.5_JPRB*PDT*(PTKEADV(JLEV) +&
       & (ZTKEPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV) )
    RHSA9(JLEV)=0.5_JPRB*PDT*(&
       & (ZAPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV) )
    RHSL9(JLEV)=0.5_JPRB*PDT*(&
       & (ZLPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV) )
    RHSI9(JLEV)=0.5_JPRB*PDT*(&
       & (ZIPT1*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV) )
  ENDDO

ELSE
  ! Don't want to extrapolate RHS, used only plain averaging
  DO JLEV=1,KFLEV
    ZU2U0(JLEV)=(PUADV(JLEV)+(PVT0(JLEV)-PVG0(JLEV))*PRCORI)*PDT
    ZV2V0(JLEV)=(PVADV(JLEV)+(PUG0(JLEV)-PUT0(JLEV))*PRCORI)*PDT
    ZT2T0(JLEV)=(PTADV(JLEV)+ PTT0(JLEV)*PKAP(JLEV)*PVVEL0(JLEV)&
     &/PRESF(JLEV))*PDT
    !ZQ2Q0(JLEV)=PDT*PQADV(JLEV)
    ZQ2Q0(JLEV)=(PQADV(JLEV)+ PQT0(JLEV)*PKAP(JLEV)*PVVEL0(JLEV)&
     &/PRESF(JLEV))*PDT
    !ZTKE2Q0(JLEV)=PDT*PTKEADV(JLEV)
    ZTKE2Q0(JLEV)=(PTKEADV(JLEV)+ PTKET0(JLEV)*PKAP(JLEV)*PVVEL0(JLEV)&
     &/PRESF(JLEV))*PDT
    ZA2A0(JLEV)=(PAT0(JLEV)*PKAP(JLEV)*PVVEL0(JLEV)&
     &/PRESF(JLEV))*PDT
    ZL2L0(JLEV)=(PLT0(JLEV)*PKAP(JLEV)*PVVEL0(JLEV)&
     &/PRESF(JLEV))*PDT
    ZI2I0(JLEV)=(PIT0(JLEV)*PKAP(JLEV)*PVVEL0(JLEV)&
     &/PRESF(JLEV))*PDT

  ! simple averaging along trajectory
  ZU2U0(JLEV)=0.5_JPRB*ZU2U0(JLEV)
  ZV2V0(JLEV)=0.5_JPRB*ZV2V0(JLEV)
  ZT2T0(JLEV)=0.5_JPRB*ZT2T0(JLEV)
  ZQ2Q0(JLEV)=0.5_JPRB*ZQ2Q0(JLEV)
  ZTKE2Q0(JLEV)=0.5_JPRB*ZTKE2Q0(JLEV)
  ZA2A0(JLEV)=0.5_JPRB*ZA2A0(JLEV)
  ZL2L0(JLEV)=0.5_JPRB*ZL2L0(JLEV)
  ZI2I0(JLEV)=0.5_JPRB*ZI2I0(JLEV)
  ! storing the same quantity for averaging
  RHSU9(JLEV)=ZU2U0(JLEV)
  RHSV9(JLEV)=ZV2V0(JLEV)
  RHST9(JLEV)=ZT2T0(JLEV)
  RHSQ9(JLEV)=ZQ2Q0(JLEV)
  RHSTKE9(JLEV)=ZTKE2Q0(JLEV)
  RHSA9(JLEV)=ZA2A0(JLEV)
  RHSL9(JLEV)=ZL2L0(JLEV)
  RHSI9(JLEV)=ZI2I0(JLEV)
  ENDDO
ENDIF



!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE ORIGIN POINT 
!              ------------------------------------

IF (NSTEP == 0) WRL9(1:KFLEV)=PETAP(1:KFLEV)

! Origin point research
CALL LARMES1C(YDGEOMETRY,YDMODEL,KFLEV,PDT,PETAP,WRL9,ZLEV)

!     ------------------------------------------------------------------

!*       3.    INTERPOLATIONS AT THE ORIGIN POINT
!              ------------------------------------


!call abor1('KONEC')
! Linear interpolation for SL physics
IF (LDSLPHY) THEN
  PSLPHY9%U(:,:)=0._JPRB
  PSLPHY9%V(:,:)=0._JPRB
  PSLPHY9%T(:,:)=0._JPRB
  PSLPHY9%GFL(:,:,:)=0._JPRB
  IF (NSTEP > 0) THEN
    ! scaling (see GPADDSLPHY)
    PSAVTEND(:,YRSLPHY%MU_SAVTEND)=PSAVTEND(:,YRSLPHY%MU_SAVTEND)*PDT
    PSAVTEND(:,YRSLPHY%MV_SAVTEND)=PSAVTEND(:,YRSLPHY%MV_SAVTEND)*PDT
    PSAVTEND(:,YRSLPHY%MT_SAVTEND)=PSAVTEND(:,YRSLPHY%MT_SAVTEND)*PDT
    PGFLSLP(:,YGFL%YQ%MPSLP)   =PGFLSLP(:,YGFL%YQ%MPSLP)*PDT
    PGFLSLP(:,YGFL%YTKE%MPSLP) =PGFLSLP(:,YGFL%YTKE%MPSLP)*PDT
    PGFLSLP(:,YGFL%YA%MPSLP)   =PGFLSLP(:,YGFL%YA%MPSLP)*PDT
    PGFLSLP(:,YGFL%YL%MPSLP)   =PGFLSLP(:,YGFL%YL%MPSLP)*PDT
    PGFLSLP(:,YGFL%YI%MPSLP)   =PGFLSLP(:,YGFL%YI%MPSLP)*PDT
!    PGFLSLP(:,YGFL%YR%MPSLP)   =PGFLSLP(:,YGFL%YR%MPSLP)*PDT 
!    PGFLSLP(:,YGFL%YS%MPSLP)   =PGFLSLP(:,YGFL%YS%MPSLP)*PDT
    CALL LARCIN1C(YDGEOMETRY,YDMODEL,KFLEV,ZLEV,PSAVTEND(:,YRSLPHY%MU_SAVTEND),PSAVTEND(:,YRSLPHY%MV_SAVTEND) &
    &  ,PSAVTEND(:,YRSLPHY%MT_SAVTEND),PGFLSLP(:,YGFL%YQ%MPSLP),PGFLSLP(:,YGFL%YTKE%MPSLP)&
    &  ,PGFLSLP(:,YGFL%YA%MPSLP),PGFLSLP(:,YGFL%YL%MPSLP),PGFLSLP(:,YGFL%YI%MPSLP) &
    &  ,PSLPHY9%U(1,:),PSLPHY9%V(1,:),PSLPHY9%T(1,:),PSLPHY9%GFL(1,:,YGFL%YQ%MP9_PH)&
    &  ,PSLPHY9%GFL(1,:,YGFL%YTKE%MP9_PH) &
    &  ,PSLPHY9%GFL(1,:,YGFL%YA%MP9_PH),PSLPHY9%GFL(1,:,YGFL%YL%MP9_PH),PSLPHY9%GFL(1,:,YGFL%YI%MP9_PH),1)
  ENDIF
ENDIF
! Linear interpolation tendency from dynamics + LS forcing 
CALL LARCIN1C(YDGEOMETRY,YDMODEL,KFLEV,ZLEV,ZQ2Q0,ZTKE2Q0,ZT2T0,ZU2U0,ZV2V0,ZA2A0,ZL2L0,ZI2I0&
             &          ,ZQL0,ZTKEL0 ,ZTL0 ,ZUL0 ,ZVL0 ,ZAL0,ZLL0,ZIL0,1)
! High order interpolation
CALL LARCIN1C(YDGEOMETRY,YDMODEL,KFLEV,ZLEV,PQT9,PTKET9,PTT9,PUT9,PVT9,PAT9,PLT9,PIT9&
             &          ,ZQT1,ZTKET1,ZTT1,ZUT1,ZVT1,ZAT1,ZLT1,ZIT1,3)

!     ------------------------------------------------------------------

!*       4.    COMPUTATION 0F TENDENCIES.
!              --------------------------

DO J=1,KFLEV
  ZTT1(J)=ZTT1(J)+ZTL0(J) +RHST9(J)
  ZQT1(J)=ZQT1(J)+ZQL0(J) +RHSQ9(J)
  ZTKET1(J)=ZTKET1(J)+ZTKEL0(J) +RHSTKE9(J)
  ZUT1(J)=ZUT1(J)+ZUL0(J) +RHSU9(J)
  ZVT1(J)=ZVT1(J)+ZVL0(J) +RHSV9(J)
  PQTEND(J)=(ZQT1(J)-PQT9(J))/PDT
  PTKETEND(J)=(ZTKET1(J)-PTKET9(J))/PDT
  PUTEND(J)=(ZUT1(J)-PUT9(J))/PDT
  PVTEND(J)=(ZVT1(J)-PVT9(J))/PDT
  PTTEND(J)=(ZTT1(J)-PTT9(J))/PDT
ENDDO

IF (LWADVCLD) THEN
  DO J=1,KFLEV
    ZAT1(J)=ZAT1(J)+ZAL0(J) +RHSA9(J)
    ZLT1(J)=ZLT1(J)+ZLL0(J) +RHSL9(J)
    ZIT1(J)=ZIT1(J)+ZIL0(J) +RHSI9(J)
    PATEND(J)=(ZAT1(J)-PAT9(J))/PDT
    PLTEND(J)=(ZLT1(J)-PLT9(J))/PDT
    PITEND(J)=(ZIT1(J)-PIT9(J))/PDT
!  WRITE(6,*) ' PLTEND',J,PLTEND(J),ZLT1(J),PLT9(J)
!  WRITE(6,*) ' PITEND',J,PITEND(J),ZIT1(J),PIT9(J)
!    PATEND(J)=0._JPRB
!    PLTEND(J)=0._JPRB
!    PITEND(J)=0._JPRB
  ENDDO

ELSE

  DO J=1,KFLEV
    PATEND(J)=0.0_JPRB
    PLTEND(J)=0.0_JPRB
    PITEND(J)=0.0_JPRB
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END SUBROUTINE LAPINE1C
