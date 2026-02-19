! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE HELDSUAREZ(YDVAB,YDCVER,YDDIMV,YDSURF, YDRIP,KIDIA, KFDIA, KLON, KLEV, &
 & PU, PV, PT, PLNSP,&
 & PGEMU,&
 & PTENU, PTENV, PTENT,&
 & PEXTRA, PEXTR2)

!**** *HELDSUAREZ * - CALCULATES SIMPLIFIED PHYSICS AND LONG TERM
!                     DIAGNOSTICS 

!     PURPOSE.
!     --------
!     - DYNAMICAL CORE EXPERIMENT
!     - THIS SUBROUTINE SUBSTITUTES THE E.C.M.W.F. PHYSICS
!       PACKAGE.

!**   Interface.
!     ----------
!        *CALL* *HELDSUAREZ*

!-----------------------------------------------------------------------

! -   INPUT ARGUMENTS.
!     -------------------

!-----------------------------------------------------------------------
! - DIMENSIONS ETC.

! KIDIA   : START OF HORIZONTAL LOOP
! KFDIA   : END OF HORIZONTAL LOOP
! KLON    : HORIZONTAL DIMENSION
! KLEV    : END OF VERTICAL LOOP AND VERTICAL DIMENSION

!-----------------------------------------------------------------------
! - FIELDS

! PU         : X-COMPONENT OF WIND.
! PV         : Y-COMPONENT OF WIND.
! PT         : TEMPERATURE.
! PLNSP      : LOG SURFACE PRESSURE

!-----------------------------------------------------------------------
! - GEOGRAPHICAL DISTRIBUTION

! PGEMU      : SINE OF LATITUDE

!-----------------------------------------------------------------------
! -   UPDATED ARGUMENTS.

! PTENU      : TENDENCY OF U-COMP. OF WIND.
! PTENV      : TENDENCY OF V-COMP. OF WIND.
! PTENT      : TENDENCY OF TEMPERATURE.
! PTENQ      : TENDENCY OF HUMIDITY
!-----------------------------------------------------------------------
!     Externals.  --- 
!     --------- 

!     Method. See documentation.
!     -------

!     Author.
!     -------
!      C. Jablonowski (ECMWF) 
!      ORIGINAL 97-10-27 Held-Suarez forcing included

!     Modifications.
!     --------------
!      08-01-08 N. Wedi, rewrite, add forcing to other state
!      K. Yessad (Jan 2011): remove useless overdimension.
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!-----------------------------------------------------------------------

USE YOMVERT            , ONLY : TVAB, VP00
USE YOMDIMV            , ONLY : TDIMV
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST             , ONLY : RPI, REA, RG, RD, RDAY, REPSM, RCPD
USE YOMLUN             , ONLY : NULOUT
USE YOMCT3             , ONLY : NSTEP
USE YOMRIP0            , ONLY : RTIMST
USE YOMRIP             , ONLY : TRIP
USE YOMCVER       , ONLY : TCVER
USE YOMDYNCORE         , ONLY : RSAMPLING_INTERVAL, RSAMPLING_START, RGAMMA_I, RGAMMA_D,&
 &                              REFERENCE_LEVEL, RPRESSURE_SCALE, L_HS_WILLIAMSON, RHS_KAPPA, LHELDSUAREZ,&
 &                              RU00_DYN, RSIGMAB, NTESTCASE, RHS_KS_KA, RPHSW_EQ, RST0_DYN, RP00_DYN,&
 &                              RT00_DYN, RDELTA_THETA, RPHSW_D_REV, RTHSW_0, RPHSW_EQ_PL_HALF, RDELTA_T,&
 &                              RHS_KA, RHS_KF, RPHSW_D
USE INTDYN_MOD         , ONLY : YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA, KLON, KLEV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNSP(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,YDSURF%YSD_XAD%NLEVS,YDSURF%YSD_XAD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,YDSURF%YSD_X2D%NDIM)
!     ------------------------------------------------------------------
!            0.2  LOCAL ARRAYS FOR THE HELD-SUAREZ TEST
!                 -------------------------------------
REAL(KIND=JPRB) :: Z_SCALED_PRESSURE
REAL(KIND=JPRB) :: Z_KV, Z_KT, Z_MAX, Z_TEQ, Z_P_I
REAL(KIND=JPRB) :: Z_SQRSIN, Z_SQRCOS
REAL(KIND=JPRB) :: Z_F1, Z_F2
INTEGER(KIND=JPIM) :: I_SAMPLING_START, I_SAMPLING_INTERVAL
INTEGER(KIND=JPIM) :: I_STAT

!--------------------------------------------------------------------
!     SETUP 
!     parameter for the Held-Suarez-Williamson forcing
!     (new radiative equilibrium temperature)
!--------------------------------------------------------------------

REAL(KIND=JPRB) :: ZEXP1, ZEXP2, ZPHI_0, ZPID180, ZHS_A
REAL(KIND=JPRB) :: ZSIGMAF(KLEV), ZSIGMAH(0:KLEV)
REAL(KIND=JPRB) :: ZSTATI, ZTIMTR, ZTETA, ZDECLIM, ZCODECM, ZNHFACT

REAL(KIND=JPRB) :: ZPRE9(KLON,0:YDDIMV%NFLEVG),ZPRE9F(KLON,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZXYB9(KLON,YDDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZST,ZN,ZN2,ZU0,ZT0,ZT02,ZVP00,ZRHO0, ZPRS, ZC, ZUU, ZUPP
REAL(KIND=JPRB) :: ZPRESHX(0:YDDIMV%NFLEVG), ZUSHEAR(YDDIMV%NFLEVG), ZTVERT(YDDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JK, JL, ITIME, IZT, IPHASE1, IPHASE2, ITEST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL         :: LLHELDSUAREZ

!     ------------------------------------------------------------------

#include "fctast.func.h"
#include "gphpre.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HELDSUAREZ',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & TSTEP=>YDRIP%TSTEP, &
 & YSD_X2D=>YDSURF%YSD_X2D, YSD_XAD=>YDSURF%YSD_XAD)
!     ------------------------------------------------------------------

IPHASE1=2400*60*60
IPHASE2=(240+24)*60*60

ZEXP1 = (RD*RGAMMA_D)/RG
ZEXP2 = (RD*RGAMMA_I)/RG
ZPID180 = RPI/180._JPRB
ZPHI_0 = 60._JPRB*ZPID180
ZHS_A = 2.65_JPRB/(15._JPRB*ZPID180)
ITIME=NINT(TSTEP)
IZT=ITIME*NSTEP
ZSTATI=REAL(IZT,JPRB)
ZTIMTR=RTIMST+ZSTATI
ZTETA=RTETA(ZTIMTR)
ZDECLIM=REL(ZTETA)
ZCODECM=COS(ZDECLIM)
ZNHFACT=(SIN(ZDECLIM)+1._JPRB)*0.5_JPRB
IF(PGEMU(KIDIA) > 0.999_JPRB) WRITE(NULOUT,*) ' ZNHFACT=',ZNHFACT

!--------------------------------------------------------------------
! compute full level pressure

DO JL=KIDIA,KFDIA
  ZPRE9 (JL,NFLEVG) = EXP(PLNSP(JL))
ENDDO
CALL GPHPRE(KLON,NFLEVG,KIDIA,KFDIA,YDVAB,YDCVER,ZPRE9,PXYB=ZXYB9,PRESF=ZPRE9F)

!=====================================================================

ZSIGMAH(0)=YDVAB%VAH(0)/VP00+YDVAB%VBH(0)
DO JK=1,KLEV
  ZSIGMAH(JK)=YDVAB%VAH(JK)/VP00+YDVAB%VBH(JK)
  ZSIGMAF(JK)=0.5_JPRB*(ZSIGMAH(JK)+ZSIGMAH(JK-1))
ENDDO

ITEST=0
LLHELDSUAREZ = LHELDSUAREZ
IF( LLHELDSUAREZ ) THEN

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA

      !=====================================================================
      !         calculate the friction for the specific layer
      !=====================================================================

      Z_SCALED_PRESSURE = ZPRE9F(JL,JK)*RPRESSURE_SCALE
      Z_MAX             = MAX (0.0_JPRB,(ZSIGMAF(JK)-RSIGMAB)/REFERENCE_LEVEL)
      Z_KV              = RHS_KF * Z_MAX
      Z_SQRSIN          = PGEMU(JL)*PGEMU(JL)
      Z_SQRCOS          = 1.0_JPRB - Z_SQRSIN

      !=====================================================================
      !         add the friction to the wind tendencies
      !=====================================================================

      PTENU(JL,JK) = PTENU(JL,JK) - Z_KV*PU(JL,JK)
      PTENV(JL,JK) = PTENV(JL,JK) - Z_KV*PV(JL,JK)

      !=====================================================================
      !         calculate the temperature relaxation, forcing function
      !=====================================================================

      Z_KT  = RHS_KA + RHS_KS_KA * Z_MAX*Z_SQRCOS*Z_SQRCOS

      !   Forcing not dependent on latitude
      !   Z_KT  = RHS_KA + RHS_KS_KA * Z_MAX

      !  IF(PGEMU(JL) > 0.0_JPRB) THEN
      !    Z_KT  = RHS_KA + RHS_KS_KA * Z_MAX*Z_SQRCOS*Z_SQRCOS
      !  ELSE
      !    Z_KT  = RHS_KA + RHS_KS_KA * Z_MAX
      !  ENDIF

      IF (L_HS_WILLIAMSON) THEN

        !=====================================================================
        !          Held-Suarez-Williamson forcing
        !=====================================================================

        IF (ZPRE9F(JL,JK) > RPHSW_D) THEN

          !!        Z_TEQ = MAX ( 200._JPRB,&
          !!         &(315._JPRB- RDELTA_T * Z_SQRSIN &
          !!         &- RDELTA_THETA * LOG(Z_SCALED_PRESSURE) * Z_SQRCOS)&
          !!         &* Z_SCALED_PRESSURE**RHS_KAPPA )
          !   Forcing not dependent on latitude in the troposphere
          !!             Z_TEQ = MAX ( 200._JPRB,&
          !!              &315._JPRB* Z_SCALED_PRESSURE**RHS_KAPPA )
          !  1 Jet in troposphere
          IF(PGEMU(JL) > 0.0_JPRB) THEN
            Z_TEQ = MAX ( 200._JPRB,&
             &(315._JPRB- RDELTA_T * Z_SQRSIN&
             &- RDELTA_THETA * LOG(Z_SCALED_PRESSURE) * Z_SQRCOS)&
             &* Z_SCALED_PRESSURE**RHS_KAPPA )
          ELSE
            Z_TEQ = MAX ( 200._JPRB,&
             &(315._JPRB&
             &- RDELTA_THETA * LOG(Z_SCALED_PRESSURE) )&
             &* Z_SCALED_PRESSURE**RHS_KAPPA )
          ENDIF

        ELSE

          !       Z_P_I = RPHSW_EQ - RPHSW_EQ_PL_HALF *&
          !        &(1.0_JPRB + TANH(ZHS_A*(ABS(ASIN(PGEMU(JL)))-ZPHI_0)))

          ! nau: HSW forcing independent on latitude in the stratosphere
          Z_P_I = RPHSW_EQ - RPHSW_EQ_PL_HALF *&
           &(1.0_JPRB + TANH(ZHS_A*(0.0_JPRB-ZPHI_0)))
          !   1 Jet in stratosphere
          !        IF(PGEMU(JL) > 0.0_JPRB) THEN
          !          Z_P_I = RPHSW_EQ - RPHSW_EQ_PL_HALF *&
          !           &(1.0_JPRB + TANH(ZHS_A*(ABS(ASIN(PGEMU(JL)))-ZPHI_0)))
          !        ELSE
          !          Z_P_I = RPHSW_EQ - RPHSW_EQ_PL_HALF *&
          !           &(1.0_JPRB + TANH(ZHS_A*(0.0_JPRB-ZPHI_0)))
          !        ENDIF
          ! nam: annual cycle on the stratospheric jets
          !        IF(PGEMU(JL) > 0.0_JPRB) THEN
          !          Z_P_I = RPHSW_EQ - RPHSW_EQ_PL_HALF *&
          !           &(1.0_JPRB + TANH(ZHS_A*(ABS(ASIN(PGEMU(JL)))*ZNHFACT-ZPHI_0)))
          !        ELSE
          !          Z_P_I = RPHSW_EQ - RPHSW_EQ_PL_HALF *&
          !           &(1.0_JPRB + TANH(ZHS_A*(ABS(ASIN(PGEMU(JL)))*(1.0_JPRB-ZNHFACT)-ZPHI_0)))
          !        ENDIF
          Z_TEQ = RTHSW_0 * (MIN(1.0_JPRB,ZPRE9F(JL,JK)*RPHSW_D_REV))**ZEXP1 &
           &+ RTHSW_0 * ( (MIN(1.0_JPRB,ZPRE9F(JL,JK)/Z_P_I))**ZEXP2 - 1.0_JPRB) 
        ENDIF

      ELSE
        !=====================================================================
        !       original Held-Suarez forcing
        !=====================================================================
        Z_TEQ = MAX ( 200._JPRB,&
         &(315._JPRB- RDELTA_T * Z_SQRSIN &
         &- RDELTA_THETA * LOG(Z_SCALED_PRESSURE) * Z_SQRCOS)&
         &* Z_SCALED_PRESSURE**RHS_KAPPA )
        !   Forcing not dependent on latitude
        !     Z_TEQ = MAX ( 200._JPRB,
        !      & 315._JPRB* Z_SCALED_PRESSURE**RHS_KAPPA )
      ENDIF
      !=====================================================================
      !         add the forcing to the temperature tendencies
      !=====================================================================

!!!!    IF( IZT <= IPHASE1) THEN

      PTENT(JL,JK) = PTENT(JL,JK) - Z_KT*(PT(JL,JK) - Z_TEQ)

    ENDDO
  ENDDO

ELSEIF( ITEST == 1 ) THEN

  !=====================================================================
  ! add forcing to attenuate to standard temperature   
  !=====================================================================
 
!  Z_KT = 1.0_JPRB/21600._JPRB
  ! upper threshold above which U=const (in meter)
  ZU0   = RU00_DYN
  ZT0   = RT00_DYN
  ZVP00 = RP00_DYN
  ZST   = RST0_DYN
  ZN    = SQRT(RG*RG/(RCPD*ZT0))
  ZRHO0 = ZVP00/(RD*ZT0)
  
  ZN2   = 2._JPRB*SQRT(RG*RG/(RCPD*ZT0))
  ZT02  = RG*RG/(RCPD*ZN2*ZN2)

  ZUPP =10500._JPRB

  ZPRESHX(0)=YDVAB%VAH(0)+YDVAB%VBH(0)*VP00
  DO JK=NFLEVG,1,-1
    ZPRESHX(JK-1)=YDVAB%VAH(JK-1)+YDVAB%VBH(JK-1)*VP00
    ZPRESHX(JK)=YDVAB%VAH(JK)+YDVAB%VBH(JK)*VP00
    ZPRS=0.5_JPRB*(ZPRESHX(JK)+ZPRESHX(JK-1))
    ZC=-(RD*ZT0/RG)*LOG(ZPRS/ZVP00)
    IF( ZC >= ZUPP ) THEN
      ZUSHEAR(JK) = ZUU
    ELSE
      ZUU = ZU0 - ZST*ZC
      ZUSHEAR(JK) = ZUU
    ENDIF
    IF( NTESTCASE == 4 ) THEN
      ZTVERT(JK) = ZT0
    ENDIF
    IF( NTESTCASE == 6 ) THEN
      IF( ZC < ZUPP ) THEN
        ZTVERT(JK) = ZT0
      ELSE
        ZTVERT(JK) = ZT02+(ZT0-ZT02)*(1.0_JPRB-TANH((ZC-ZUPP)/3000._JPRB));
      ENDIF
    ENDIF

  ENDDO

  ! attenuate to zonal isothermal flow with shear 
  Z_KT = 1.0_JPRB/100._JPRB
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PTENU(JL,JK) = PTENU(JL,JK) - Z_KT*(PU(JL,JK) - ZUSHEAR(JK))
      PTENV(JL,JK) = PTENV(JL,JK) - Z_KT*PV(JL,JK)
      PTENT(JL,JK) = PTENT(JL,JK) - Z_KT*(PT(JL,JK) - ZTVERT(JK))
    ENDDO
  ENDDO

ELSEIF( ITEST == 2 ) THEN

  !============================================================================
  ! add vertical Laplacian forcing (sea breeze) with prescribed boundary values
  !============================================================================
 
  ZU0   = RU00_DYN
  ZT0   = RT00_DYN
  ZVP00 = RP00_DYN
  ZST   = RST0_DYN
  ZN    = SQRT(RG*RG/(RCPD*ZT0))
  ZRHO0 = ZVP00/(RD*ZT0)

  ! CALL SPTSV
  !if (JPRB == N_DEFAULT_REAL_KIND) then
  !  call SPTSV (iter,1,zwork(1,1),zwork(2,2),zwork(1,3),kmaxit+1,info)
  !ELSEIF (JPRB == N_DOUBLE_KIND) then
  !  call DPTSV (iter,1,zwork(1,1),zwork(2,2),zwork(1,3),kmaxit+1,info)
  !else
  !  call ABOR1 ('REAL(KIND=JPRB) is neither default real nor double precision')
  !endif 
  !DO JK=1,KLEV
  !  DO JL=KIDIA,KFDIA
  !    PTENT(JL,JK) = PTENT(JL,JK) + zwork 
  !  ENDDO
  !ENDDO

  !=====================================================================
  !     approximate normal mode disturbance
  !=====================================================================
  
!!!!    IF( IPHASE1 < IZT.AND.IZT <= IPHASE2) THEN
  !!      PEXTRA (JL,JK,2) =  Z_U_D 
  !
  ! add the disturbance to the wind tendencies
  !
  !!      IF(PGEMU(JL) > 0.0_JPRB) THEN    
  !!        PTENU(JL,JK) = PTENU(JL,JK) - Z_U_D
  !!      ENDIF
  !!    ELSE
  !!      PEXTRA (JL,JK,2) = 0._JPRB
!!!!    ENDIF

ENDIF
  
!=======================================================================
!     DIAGNOSTICS:: sampling period is day SAMPLING_START - open end
!=======================================================================

!     put the initialization of the two sampling variables
!     I_SAMPLING_START and I_SAMPLING_INTERVAL as well as
!     the initialization of the extra fields in 
!     a SETUP routine later on !!!! 

I_SAMPLING_START    = NINT(RSAMPLING_START/TSTEP)
I_SAMPLING_INTERVAL = NINT(RSAMPLING_INTERVAL/TSTEP)

!!IF ( NSTEP == 0 ) THEN
!!  DO JK=1,KLEV
!!    DO JL=KIDIA,KFDIA
!!      PEXTRA (JL,JK,1) = 0.0_JPRB
!!      PEXTRA (JL,JK,2) = 0.0_JPRB
!!      PEXTRA (JL,JK,3) = 0.0_JPRB
!!      PEXTRA (JL,JK,4) = 0.0_JPRB
!!      PEXTRA (JL,JK,5) = 0.0_JPRB
!!      PEXTRA (JL,JK,6) = 0.0_JPRB
!!      PEXTRA (JL,JK,7) = 0.0_JPRB
!!      PEXTRA (JL,JK,8) = 0.0_JPRB
!!      PEXTRA (JL,JK,9) = 0.0_JPRB
!!      PEXTRA (JL,JK,10) = 0.0_JPRB
!!      PEXTRA (JL,JK,11) = 0.0_JPRB
!!      PEXTRA (JL,JK,12) = 0.0_JPRB
!!    ENDDO
!!  ENDDO
!!  DO JL=KIDIA,KFDIA
!!    PEXTR2 (JL,1) = 0.0_JPRB
!!    PEXTR2 (JL,2) = 0.0_JPRB
!!  ENDDO

!!  WRITE (NULOUT,*) ' SUBROUTINE HELDSUAREZ, INITIALIZATION'

!!ENDIF

!     sampling interval day SAMPLING_START - open end          

IF (     NSTEP >= I_SAMPLING_START &
   &.AND. MOD((NSTEP-I_SAMPLING_START),I_SAMPLING_INTERVAL) == 0)&
&THEN

!=======================================================================
!     calculate the time mean values
!=======================================================================

! add special diagnostics here

!=======================================================================
!       coefficients to calculate the time mean over i_stat timesteps
!=======================================================================

  I_STAT = (NSTEP-I_SAMPLING_START)/I_SAMPLING_INTERVAL + 1
  Z_F1 = 1.0_JPRB/REAL(I_STAT,JPRB)
  Z_F2 = REAL(I_STAT-1,JPRB)*Z_F1

!=======================================================================
!       multi level fields
!=======================================================================

!!  DO JK=1,KLEV
!!    DO JL=KIDIA,KFDIA
!!      PEXTRA (JL,JK,1) =Z_F1*PU(JL,JK)+Z_F2*PEXTRA(JL,JK,1)
!!      PEXTRA (JL,JK,2) =Z_F1*PV(JL,JK)+Z_F2*PEXTRA(JL,JK,2)
!!      PEXTRA (JL,JK,3) =Z_F1*PT(JL,JK)+Z_F2*PEXTRA(JL,JK,3)
!!      PEXTRA (JL,JK,4) =Z_F1*PAPHIF(JL,JK)+Z_F2*PEXTRA(JL,JK,4)
!!      PEXTRA (JL,JK,5) =Z_F1*PVERVEL(JL,JK)+Z_F2*PEXTRA(JL,JK,5)
!!      PEXTRA (JL,JK,6) =Z_F1*PU(JL,JK)*PU(JL,JK) +Z_F2*PEXTRA(JL,JK,6)
!!      PEXTRA (JL,JK,7) =Z_F1*PU(JL,JK)*PV(JL,JK) +Z_F2*PEXTRA(JL,JK,7)
!!      PEXTRA (JL,JK,8) =Z_F1*PV(JL,JK)*PT(JL,JK) +Z_F2*PEXTRA(JL,JK,8)
!!      PEXTRA (JL,JK,9) =Z_F1*PV(JL,JK)*PV(JL,JK) +Z_F2*PEXTRA(JL,JK,9)
!!      PEXTRA (JL,JK,10) =Z_F1*PT(JL,JK)*PT(JL,JK) +Z_F2*PEXTRA(JL,JK,10)
!!      PEXTRA (JL,JK,11) =Z_F1*PVERVEL(JL,JK)*&
!!       &PVERVEL(JL,JK) +Z_F2*PEXTRA(JL,JK,11)
!!      PEXTRA (JL,JK,12) =Z_F1*PVERVEL(JL,JK)*&
!!       &PT(JL,JK) +Z_F2*PEXTRA(JL,JK,12)
!!    ENDDO
!!  ENDDO

!=======================================================================
!       single level field ps
!=======================================================================

!       2-d extra fields

!!  DO JL=KIDIA,KFDIA
!!    PEXTR2 (JL,1) =Z_F1*PAPRS (JL,KLEV)+Z_F2*PEXTR2(JL,1)
!!    PEXTR2 (JL,2) =Z_F1*PAPRS (JL,KLEV)*PAPRS (JL,KLEV)+Z_F2*PEXTR2(JL,2)
!!  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('HELDSUAREZ',1,ZHOOK_HANDLE)
END SUBROUTINE HELDSUAREZ
