! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE CLOUD_SATADJ(YDECLDP, YDEPHLI, YDECUMF, YDEPHY, YDSPP_CONFIG,&
 & KIDIA,    KFDIA,    KLON,    KLEV,  KTYPE, &
 & PTSPHY, PAP,  PAPH, &
 & PT, PQ, PA, &
 & PL, PI, & 
 & PTENDENCY_CML_T, PTENDENCY_CML_Q, PTENDENCY_CML_A, &
 & PTENDENCY_CML_L, PTENDENCY_CML_I, &
 & PHRSW,    PHRLW, &
 & PTENDENCY_VDF_T, PTENDENCY_VDF_Q, &
 & PMFU,     PMFD, PLUDELI, &
 & PVERVEL,  PGP2DSPP, &
 & PTENDENCY_LOC_T, PTENDENCY_LOC_Q, PTENDENCY_LOC_A, &
 & PTENDENCY_LOC_L, PTENDENCY_LOC_I, &
 & PEXTRA,   KFLDX)
 
!===============================================================================
!**** *CLOUD_SATADJ* -  ROUTINE FOR PARAMATERIZATION OF CLOUD PROCESSES
!                  FOR PROGNOSTIC CLOUD SCHEME
!!
!     R.Forbes     (E.C.M.W.F.)
!!
!     PURPOSE
!     -------
!     Evaporation/condensation of cloud water by adjustment to saturation
!     changed by diabatic heating/cooling from ascent/subsidence/radiation
!!
!     INTERFACE.
!     ----------
!     *CLOUD_SATADJ* is called from *CALLPAR*
!     Input tendencies from dynamics (omega), radiation, turbulent mixing,
!     and the convection scheme. Outputs changes to temperature, humidity, 
!     cloud fraction and cloud condensate.
!!
!     EXTERNALS.
!     ----------
!          QSATWATADJ, QSATMIXADJ
!!
!     MODIFICATIONS.
!     -------------
!     01-10-2016 : R.Forbes  New routine. Duplicate of cond/evap from cloudsc.F90
!     15-02-2020 : R.Forbes  Rewrite of routine, modified saturation adjustment and limits 
!                            New simpler fn QSATMIXADJ/QSATWATADJ to replace CUADJTQ
!
!     REFERENCES.
!     ----------
!     Tietdke MWR 1993 - original description of the cloud parametrization
!     Jakob PhD 2000
!     Gregory et al. (2000) QJRMS
!     Tompkins el al. (2007) QJRMS - ice supersaturation parametrization
!!
!===============================================================================

USE YOECLDP  , ONLY : TECLDP
USE YOEPHLI  , ONLY : TEPHLI
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RTT, RLVTT, RLSTT
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2
USE YOECUMF  , ONLY : TECUMF
USE YOEPHY   , ONLY : TEPHY
USE SPP_MOD     , ONLY : TSPP_CONFIG
USE SPP_GEN_MOD , ONLY : SPP_PERT

IMPLICIT NONE

!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------
 
TYPE(TECLDP)      ,INTENT(IN)    :: YDECLDP
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
TYPE(TSPP_CONFIG) ,INTENT(IN)    :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON)      ! Convection type 0-3
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY           ! Physics timestep
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)   ! Pressure on full levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
! Initial state
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)    ! Temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)    ! Humidity (kg/kg)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV)    ! Cloud fraction (0-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV)    ! Cloud liquid (kg/kg)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV)    ! Cloud ice (kg/kg)
! Input tendency to add to initial state
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_T(KLON,KLEV) ! Temperature
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_Q(KLON,KLEV) ! Humidity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_A(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_L(KLON,KLEV) ! Cloud liquid
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_CML_I(KLON,KLEV) ! Cloud ice
! Input tendencies from other processes
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV)   ! Short-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV)   ! Long-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_VDF_T(KLON,KLEV) ! Vert diff T
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDENCY_VDF_Q(KLON,KLEV) ! Vert diff Q
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)    ! Conv. mass flux up
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV)    ! Conv. mass flux down
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLUDELI(KLON,KLEV,4)! Conv. detrained liq/ice/vapor/T
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) ! Vertical velocity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YDSPP_CONFIG%SM%NRFTOTAL) ! perturbation pattern
! Local output tendencies from cloud scheme
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_T(KLON,KLEV) ! Temperature
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_Q(KLON,KLEV) ! Humidity
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_A(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_L(KLON,KLEV) ! Cloud liquid
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDENCY_LOC_I(KLON,KLEV) ! Cloud ice
! Extra fields for diagnostics
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX ! Number of extra fields

!-------------------------------------------------------------------------------
!                       Declare local variables
!-------------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZT(KLON)   
REAL(KIND=JPRB) :: ZQ(KLON)
REAL(KIND=JPRB) :: ZA(KLON)
REAL(KIND=JPRB) :: ZL(KLON)
REAL(KIND=JPRB) :: ZI(KLON)

REAL(KIND=JPRB) :: ZLI(KLON)
REAL(KIND=JPRB) :: ZT_UPD(KLON)   
REAL(KIND=JPRB) :: ZQ_UPD(KLON)   
REAL(KIND=JPRB) :: ZA_UPD(KLON)   

REAL(KIND=JPRB) :: ZLIQFRAC(KLON) ! cloud liquid water fraction: ql/(ql+qi)
REAL(KIND=JPRB) :: ZICEFRAC(KLON) ! cloud ice water fraction: qi/(ql+qi)
REAL(KIND=JPRB) :: ZQSMIX(KLON)   ! diagnostic mixed phase saturation 
REAL(KIND=JPRB) :: ZQSLIQ(KLON)   ! liquid water saturation
REAL(KIND=JPRB) :: ZQSICE(KLON)   ! ice water saturation
REAL(KIND=JPRB) :: ZQSLIQK(KLON)  ! water saturation to T=-40, Koop T<-40
REAL(KIND=JPRB) :: ZQSLIM(KLON)
REAL(KIND=JPRB) :: ZFOEEICE(KLON)
REAL(KIND=JPRB) :: ZFOEELIQ(KLON)
REAL(KIND=JPRB) :: ZFOEEMIX(KLON)
REAL(KIND=JPRB) :: ZFOEALFA(KLON)
!REAL(KIND=JPRB) :: ZADVW(KLON), ZADVWD(KLON)

REAL(KIND=JPRB) :: ZACOND(KLON)
REAL(KIND=JPRB) :: ZLCOND_CLD(KLON) 
REAL(KIND=JPRB) :: ZLCOND_ENV(KLON)
REAL(KIND=JPRB) :: ZLCOND_CLD_L(KLON) 
REAL(KIND=JPRB) :: ZLCOND_ENV_L(KLON)
REAL(KIND=JPRB) :: ZLCOND_CLD_I(KLON) 
REAL(KIND=JPRB) :: ZLCOND_ENV_I(KLON)
REAL(KIND=JPRB) :: ZLEVAP_CLD_L(KLON)
REAL(KIND=JPRB) :: ZLEVAP_CLD_I(KLON)
REAL(KIND=JPRB) :: ZLEVAP_CLD(KLON)  
REAL(KIND=JPRB) :: ZAEVAP(KLON)
REAL(KIND=JPRB) :: ZFOKOOP(KLON)
REAL(KIND=JPRB) :: ZLICLD(KLON)
REAL(KIND=JPRB) :: ZLIQCLD(KLON)
REAL(KIND=JPRB) :: ZICECLD(KLON)
REAL(KIND=JPRB) :: ZDQS(KLON)
REAL(KIND=JPRB) :: ZTTMP(KLON)
REAL(KIND=JPRB) :: ZQTMP(KLON)  
REAL(KIND=JPRB) :: ZDTGDP(KLON) 
REAL(KIND=JPRB) :: ZRDTGDP(KLON)
REAL(KIND=JPRB) :: ZRHO(KLON)
REAL(KIND=JPRB) :: ZGDP(KLON)
REAL(KIND=JPRB) :: ZDA(KLON)
REAL(KIND=JPRB) :: ZDP(KLON)
REAL(KIND=JPRB) :: ZDZ(KLON)
REAL(KIND=JPRB) :: ZXRAMID
REAL(KIND=JPRB) :: ZP_R(KLON)
REAL(KIND=JPRB) :: ZSUPSATL(KLON)
REAL(KIND=JPRB) :: ZSUPSATI(KLON)
REAL(KIND=JPRB) :: ZSUPSATA(KLON)
REAL(KIND=JPRB) :: ZCORQSLIQ(KLON)
REAL(KIND=JPRB) :: ZCORQSMIX(KLON)
REAL(KIND=JPRB) :: ZEVAPLIMMIX(KLON)
REAL(KIND=JPRB) :: ZQSMIX_UPD(KLON)   ! diagnostic mixed phase saturation 
REAL(KIND=JPRB) :: ZQSLIQ_UPD(KLON)   ! liquid water saturation
REAL(KIND=JPRB) :: ZQSICE_UPD(KLON)   ! ice water saturation
REAL(KIND=JPRB) :: ZQSLIQK_UPD(KLON)  ! water saturation to T=-40, Koop T<-40
REAL(KIND=JPRB) :: ZCORQSMIX_UPD(KLON)
REAL(KIND=JPRB) :: ZPVERVEL(KLON)
REAL(KIND=JPRB) :: ZDTFORCE(KLON)
REAL(KIND=JPRB) :: ZDQFORCE(KLON)
REAL(KIND=JPRB) :: ZT_ADJ(KLON) ! Supersat check temperature change
REAL(KIND=JPRB) :: ZQ_ADJ(KLON) ! Supersat check humidity change
REAL(KIND=JPRB) :: ZA_ADJ(KLON) ! Supersat check cloud fraction change
REAL(KIND=JPRB) :: ZL_ADJ(KLON) ! Supersat check cloud water change
REAL(KIND=JPRB) :: ZI_ADJ(KLON) ! Supersat check cloud ice change

REAL(KIND=JPRB) :: ZQSATNEW
REAL(KIND=JPRB) :: ZANEW
REAL(KIND=JPRB) :: ZDTDIAB
!REAL(KIND=JPRB) :: ZMFDN
REAL(KIND=JPRB) :: ZCOR
REAL(KIND=JPRB) :: ZCDMAX
REAL(KIND=JPRB) :: ZLCONDLIM
REAL(KIND=JPRB) :: ZDPMXDT
REAL(KIND=JPRB) :: ZDTDP
REAL(KIND=JPRB) :: ZEPSEC
REAL(KIND=JPRB) :: ZFAC
REAL(KIND=JPRB) :: ZFACW
REAL(KIND=JPRB) :: ZQE
REAL(KIND=JPRB) :: ZQTMST
REAL(KIND=JPRB) :: ZRDCP
REAL(KIND=JPRB) :: ZRHC(KLON)
REAL(KIND=JPRB) :: ZSIGK
REAL(KIND=JPRB) :: ZWTOT
REAL(KIND=JPRB) :: ZZDL
REAL(KIND=JPRB) :: ZZZDT
REAL(KIND=JPRB) :: ZTMPA
REAL(KIND=JPRB) :: ZEPSILON

! A bunch of SPP variables 
LOGICAL            :: LLPERT_QSATVERVEL, LLPERT_RAMID
INTEGER(KIND=JPIM) :: IPRAMID, IPQSATVERVEL
INTEGER(KIND=JPIM) :: IPN ! SPP perturbation pointer
TYPE(SPP_PERT)     :: PN1QSATVERVEL, PN1RAMID

! Controls T-dependence for liquid/ice production (1=mixedphase, 2=homogfrz)
INTEGER(KIND=JPIM) :: IFTLIQICE 

INTEGER(KIND=JPIM) :: JL, JK, IK, IS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cloud_supersatcheck.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "fccld.func.h"

!===============================================================================
IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ',0,ZHOOK_HANDLE)

ASSOCIATE(NCLDTOP=>YDECLDP%NCLDTOP, &
 & NSSOPT=>YDECLDP%NSSOPT, & 
 & LBUD23=>YDEPHY%LBUD23, &
 & LCLDBUDC=>YDECLDP%LCLDBUDC, &
 & LCLDBUDL=>YDECLDP%LCLDBUDL, &
 & LCLDBUDI=>YDECLDP%LCLDBUDI, &
 & LCLDBUDT=>YDECLDP%LCLDBUDT, &
 & LCLDBUD_VERTINT=>YDECLDP%LCLDBUD_VERTINT, &
 & LCLDBUD_TIMEINT=>YDECLDP%LCLDBUD_TIMEINT, &
 & RAMID=>YDECLDP%RAMID,   & 
 & RAMIN=>YDECLDP%RAMIN,   &
 & RLMIN=>YDECLDP%RLMIN,   &
 & RTHOMO=>YDECLDP%RTHOMO, &
 & RKOOPTAU=>YDECLDP%RKOOPTAU, &
 & RMFADVW=>YDECUMF%RMFADVW, RMFADVWDD=>YDECUMF%RMFADVWDD )

!===============================================================================


!######################################################################
!
!             1.  *** INITIAL VALUES FOR VARIABLES ***
!
!######################################################################

! Define a small number
ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

! ---------------------
! Some simple constants
! ---------------------
ZQTMST  = 1.0_JPRB/PTSPHY
ZRDCP   = RD/RCPD
ZEPSEC  = 1.E-14_JPRB

! -----------------------------------------------------------
! Parameters for convective subsidence treatment in dynamics
! (currently commented out as not active in this routine)
! -----------------------------------------------------------
!DO JL=KIDIA, KFDIA
!  ZADVW(JL)=1.0_JPRB
!  ZADVWD(JL)=1.0_JPRB
!  IF(KTYPE(JL)==1.AND.RMFADVW>0) THEN
!    ZADVW(JL) =1.0_JPRB-RMFADVW
!    ZADVWD(JL)=1.0_JPRB-RMFADVW*RMFADVWDD
!  ENDIF
!ENDDO

! ------------------------------------------------------------------
! SPP perturbations to adiabatic vertical velocity that affects saturation adjustment
! ------------------------------------------------------------------
IF (YDSPP_CONFIG%LSPP) THEN
  ! Saturation due to adiabatic vertical velocity
  !
  IPN = YDSPP_CONFIG%PPTR%QSATVERVEL
  LLPERT_QSATVERVEL= IPN > 0
  IF (LLPERT_QSATVERVEL) THEN
    PN1QSATVERVEL  = YDSPP_CONFIG%SM%PN(IPN)
    IPQSATVERVEL   = PN1QSATVERVEL%MP
  ENDIF

  ! Critical relative humidity
  !
  IPN = YDSPP_CONFIG%PPTR%RAMID
  LLPERT_RAMID= IPN > 0
  IF (LLPERT_RAMID) THEN
    PN1RAMID  = YDSPP_CONFIG%SM%PN(IPN)
    IPRAMID   = PN1RAMID%MP
  ENDIF
ELSE
  LLPERT_RAMID      =.FALSE.
  LLPERT_QSATVERVEL =.FALSE.
ENDIF


! ------------------------------------------------------------------
! Zero tendency arrays at model top levels where no cloud processes 
! ------------------------------------------------------------------
DO JK=1,NCLDTOP-1
  DO JL=KIDIA,KFDIA
    PTENDENCY_LOC_T(JL,JK) = 0.0_JPRB
    PTENDENCY_LOC_Q(JL,JK) = 0.0_JPRB 
    PTENDENCY_LOC_A(JL,JK) = 0.0_JPRB
    PTENDENCY_LOC_L(JL,JK) = 0.0_JPRB
    PTENDENCY_LOC_I(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!----------------------------------------------------------------------
!
!                   START OF VERTICAL LOOP OVER LEVELS
!
!----------------------------------------------------------------------
DO JK=NCLDTOP,KLEV
  
  ! ------------------------------
  ! Updated state initialization 
  ! ------------------------------
  DO JL=KIDIA,KFDIA
    ZT(JL)  = PT(JL,JK) + PTSPHY*PTENDENCY_CML_T(JL,JK)
    ZQ(JL)  = PQ(JL,JK) + PTSPHY*PTENDENCY_CML_Q(JL,JK) 
    ZA(JL)  = PA(JL,JK) + PTSPHY*PTENDENCY_CML_A(JL,JK)
    ZL(JL)  = PL(JL,JK) + PTSPHY*PTENDENCY_CML_L(JL,JK)
    ZI(JL)  = PI(JL,JK) + PTSPHY*PTENDENCY_CML_I(JL,JK)
    ! Ensure cloud fraction is between 0 and 1
    ZA(JL)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZA(JL)))
  ENDDO
  
  ! ------------------------------------------------------------------
  ! Perturb vertical velocity for adiabatic ascent forcing for condensation 
  ! ------------------------------------------------------------------
  DO JL=KIDIA,KFDIA
    IF (LLPERT_QSATVERVEL) THEN  ! Compute SPP perturbation:
      ZPVERVEL(JL)=PVERVEL(JL,JK)*EXP(PN1QSATVERVEL%MU(1)+ &
              & PN1QSATVERVEL%XMAG(1)*PGP2DSPP(JL,IPQSATVERVEL))
    ELSE
      ZPVERVEL(JL)=PVERVEL(JL,JK)  ! (unperturbed)
    ENDIF
  ENDDO

  ! ------------------------------------------------------------------
  ! Define T,Q forcing terms for condensation/evaporation 
  ! ------------------------------------------------------------------
  DO JL=KIDIA,KFDIA

    ZDP(JL)     = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
    ZGDP(JL)    = RG/ZDP(JL)                    ! g/dp
    ZDTGDP(JL)  = PTSPHY*ZGDP(JL)               ! dt g/dp
    ZRDTGDP(JL) = ZDP(JL)*(1.0_JPRB/(PTSPHY*RG))! 1/(dt g/dp)

    ZDTDP   = ZRDCP*PT(JL,JK)/PAP(JL,JK)
    ZDPMXDT = ZDP(JL)*ZQTMST

    ! Vertical velocity from dynamics (omega) and convective subsidence
    ! (latter commented out due to double counting of subsid evap in cloudsc) 
    !ZMFDN   = 0.0_JPRB
    !IF(JK < KLEV) ZMFDN=PMFU(JL,JK+1)+PMFD(JL,JK+1)
    !ZWTOT   = ZPVERVEL(JL)+(0.5_JPRB*RG*(PMFU(JL,JK)*ZADVW(JL)+PMFD(JL,JK)*ZADVWD(JL)+ZMFDN))
    
    ZWTOT   = ZPVERVEL(JL)
    ZWTOT   = MIN(ZDPMXDT,MAX(-ZDPMXDT,ZWTOT))

    ! SW and LW radiative heating/cooling
    ZZZDT   = PHRSW(JL,JK)+PHRLW(JL,JK)
    ZDTDIAB = MIN(ZDPMXDT*ZDTDP,MAX(-ZDPMXDT*ZDTDP,ZZZDT))*PTSPHY

    ! Temperature forcing term
    ZDTFORCE(JL) = ZDTDP*ZWTOT*PTSPHY + ZDTDIAB + PLUDELI(JL,JK,4)*ZDTGDP(JL) &
                                    & + PTENDENCY_VDF_T(JL,JK)*PTSPHY

    ! Humidity forcing term
    ZDQFORCE(JL) = PLUDELI(JL,JK,3)*ZDTGDP(JL) + PTENDENCY_VDF_Q(JL,JK)*PTSPHY

  ENDDO

  DO JL=KIDIA,KFDIA

    !----------------------------------------------------------------------
    !
    ! 1. INITIALIZE VARIABLES
    !
    !----------------------------------------------------------------------

    !-------------------------
    ! derived variables needed
    !-------------------------

    ZDP(JL)  = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
    ZRHO(JL) = PAP(JL,JK)/(RD*PT(JL,JK))     ! p/RT air density
    ZDZ(JL)  = ZDP(JL)/(ZRHO(JL)*RG)         ! layer depth (m)
    ZP_R(JL) = 1.0_JPRB/PAP(JL,JK) 

    !-------------------------------------------------------------------
    ! Calculate liq/ice fractions (no longer a diagnostic relationship)
    !-------------------------------------------------------------------
    ZLI(JL)=ZL(JL)+ZI(JL)
    IF (ZLI(JL)>RLMIN) THEN
      ZLIQFRAC(JL)=ZL(JL)/ZLI(JL)
      ZICEFRAC(JL)=1.0_JPRB-ZLIQFRAC(JL)
    ELSE
      ZLIQFRAC(JL)=0.0_JPRB
      ZICEFRAC(JL)=0.0_JPRB
    ENDIF

    !----------------------------------------------
    ! in-cloud consensate amount for updated state
    !----------------------------------------------
    ZTMPA       = 1.0_JPRB/MAX(ZA(JL),0.01_JPRB)
    ZLIQCLD(JL) = ZL(JL)*ZTMPA
    ZICECLD(JL) = ZI(JL)*ZTMPA
    ZLICLD(JL)  = ZLIQCLD(JL)+ZICECLD(JL)

    !---------------------------
    ! Critical relative humidity
    !---------------------------
    IF (LLPERT_RAMID) THEN !Apply SPP perturbations
      ZXRAMID= RAMID*EXP(PN1RAMID%MU(1)+PN1RAMID%XMAG(1)*PGP2DSPP(JL, IPRAMID))
      IF (YDSPP_CONFIG%LRAMIDLIMIT1) THEN
        ZXRAMID= MIN(1.0_JPRB, ZXRAMID)
      ENDIF
    ELSE
      ZXRAMID= RAMID ! (unperturbed)
    ENDIF

    ZRHC(JL)=ZXRAMID
    ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
    ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
    IF(ZSIGK > 0.8_JPRB) THEN
      ZRHC(JL)=ZXRAMID+(1.0_JPRB-ZXRAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
    ENDIF

    !--------------------------------------------------
    ! Saturation functions for start of timestep state
    !--------------------------------------------------

    ! Saturation wrt water (all T)
    ZFOEELIQ(JL) = MIN(FOEELIQ(PT(JL,JK))*ZP_R(JL),0.5_JPRB)
    ZQSLIQ(JL)   = ZFOEELIQ(JL)/(1.0_JPRB-RETV*ZFOEELIQ(JL))

    ! Liquid-phase saturation adjustment factor to bring gridbox T,Q to saturation
    ZFACW       = R5LES/((PT(JL,JK)-R4LES)**2)
    ZCOR        = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEELIQ(JL))
    ZCORQSLIQ(JL)   = 1.0_JPRB+RALVDCP*ZFACW*ZCOR*ZQSLIQ(JL)

    ! Saturation wrt ice (T<0), water (T>0)
    ZFOEEICE(JL) = MIN(FOEEW(PT(JL,JK))*ZP_R(JL),0.5_JPRB)
    ZQSICE(JL)   = ZFOEEICE(JL)/(1.0_JPRB-RETV*ZFOEEICE(JL))

    ! Saturation wrt Koop curve (T<0) liquid water saturation limited to Koop curve
    ! FOEDELTA = 1 for T>0, =0 for T<0
    ZFOKOOP(JL)  = MIN(ZFOEEICE(JL)*(RKOOP1-RKOOP2*PT(JL,JK)),ZFOEELIQ(JL))
    ZQSLIQK(JL)  = ZFOKOOP(JL)/(1.0_JPRB-RETV*ZFOKOOP(JL))

    ! Saturation wrt water (T>0), mixed (RTT<T<0), ice (T<RTT)
    ZFOEALFA(JL) = FOEALFA(PT(JL,JK))
    ZFOEEMIX(JL) = MIN(FOEEWM(PT(JL,JK))*ZP_R(JL),0.5_JPRB)
    ZCOR         = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEMIX(JL))
    ZQSMIX(JL)   = ZFOEEMIX(JL)*ZCOR

    ! Reduce maximum supersaturation to between water/Koop saturation and ice saturation
    ! QSLIQK = QSMIX = QSLIQ for T>0
    ZQSLIQK(JL) = 0.5_JPRB*(ZQSLIQK(JL)+ZQSMIX(JL))

    ! Mixed-phase saturation adjustment factor to bring gridbox T,Q to saturation
    ZCORQSMIX(JL) = 1.0_JPRB+ZQSMIX(JL)*ZCOR*FOEDEM(PT(JL,JK))

    ! Evaporation/sublimation limits
    ZEVAPLIMMIX(JL) = MAX((ZQSMIX(JL)-PQ(JL,JK))/ZCORQSMIX(JL),0.0_JPRB)

    ! Calculate the maximum grid-box mean saturated state taking into account
    ! the saturation assumption in-cloud and ice supersaturated Koop limit
    ! in the clear air. Note ZQSLIQK = ZQSMIX for T warmer than 0degC
    ZQSLIM(JL) = PA(JL,JK)*ZQSMIX(JL) + (1.0_JPRB-PA(JL,JK))*ZQSLIQK(JL)

    !--------------------------------------------------
    ! Saturation functions for updated state
    !--------------------------------------------------

    ! Saturation wrt water (all T)
    ZFOEELIQ(JL)    = MIN(FOEELIQ(ZT(JL))*ZP_R(JL),0.5_JPRB)
    ZQSLIQ_UPD(JL)  = ZFOEELIQ(JL)/(1.0_JPRB-RETV*ZFOEELIQ(JL))

    ! Saturation wrt ice (T<0), water (T>0)
    ZFOEEICE(JL)    = MIN(FOEEW(ZT(JL))*ZP_R(JL),0.5_JPRB)
    ZQSICE_UPD(JL)  = ZFOEEICE(JL)/(1.0_JPRB-RETV*ZFOEEICE(JL))

    ! Saturation wrt Koop curve (T<0) liquid water saturation limited to Koop curve
    ! FOEDELTA = 1 for T>0, =0 for T<0
    ZFOKOOP(JL)     = MIN(ZFOEEICE(JL)*(RKOOP1-RKOOP2*ZT(JL)),ZFOEELIQ(JL))
    ZQSLIQK_UPD(JL) = ZFOKOOP(JL)/(1.0_JPRB-RETV*ZFOKOOP(JL))

    ! Saturation wrt water (T>0), mixed (RTT<T<0), ice (T<RTT)
    ZFOEALFA(JL)    = FOEALFA(ZT(JL))
    ZFOEEMIX(JL)    = MIN(FOEEWM(ZT(JL))*ZP_R(JL),0.5_JPRB)
    ZCOR            = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEMIX(JL))
    ZQSMIX_UPD(JL)  = ZFOEEMIX(JL)*ZCOR

    ! Reduce maximum supersaturation to between water/Koop saturation and ice saturation
    ! QSLIQK = QSMIX = QSLIQ for T>0
    ZQSLIQK_UPD(JL) = 0.5_JPRB*(ZQSLIQK_UPD(JL)+ZQSMIX_UPD(JL))

    ! Mixed-phase saturation adjustment factor to bring gridbox T,Q to saturation
    ZCORQSMIX_UPD(JL) = 1.0_JPRB+ZQSMIX_UPD(JL)*ZCOR*FOEDEM(ZT(JL))

    ! Evaporation/sublimation limits
    !ZEVAPLIMMIX_UPD(JL) = MAX((ZQSMIX_UPD(JL)-ZQ(JL))/ZCORQSMIX_UPD(JL),0.0_JPRB)

    ! Calculate the maximum grid-box mean saturated state taking into account
    ! the saturation assumption in-cloud and ice supersaturated Koop limit
    ! in the clear air. Note ZQSLIQK = ZQSMIX for T warmer than 0degC
    !ZQSLIM_UPD(JL) = ZA(JL)*ZQSMIX_UPD(JL) + (1.0_JPRB-ZA(JL))*ZQSLIQK_UPD(JL)

  ENDDO
    
  !======================================================================
  !
  !
  !  2.  CONDENSATION/EVAPORATION DUE TO DQSAT/DT
  !
  !
  !======================================================================
  !  calculate dqs/dt
  !  Note: For the separate prognostic Qi and Ql, one would ideally use
  !  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
  !  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
  !  These would then instantaneous freeze if T<-38C or lead to ice growth 
  !  by deposition in warmer mixed phase clouds.  However, since we do 
  !  not have a separate prognostic equation for in-cloud humidity or a 
  !  statistical scheme approach in place, the depositional growth of ice 
  !  in the mixed phase can not be modelled and we resort to supersaturation  
  !  wrt ice instanteously converting to ice over one timestep 
  !  (see Tompkins et al. QJRMS 2007 for details)
  !  Thus for the initial implementation the diagnostic mixed phase is 
  !  retained for the moment, and the level of approximation noted.  
  !----------------------------------------------------------------------

  DO JL=KIDIA,KFDIA

    ! Add forcing term to start of timestep temperature
    ZTTMP(JL) = PT(JL,JK) + ZDTFORCE(JL)
    ZTTMP(JL) = MAX(ZTTMP(JL),160.0_JPRB) ! security

    ! Add forcing term to start of timestep saturation specific humidity
    ZQTMP(JL) = ZQSMIX(JL) + ZDQFORCE(JL)

    ! Perform saturation adjustment
    ZQSATNEW = QSATMIXADJ(ZTTMP(JL), ZQTMP(JL), PAP(JL,JK))

    ! Calculate the condensation/evaporation as the change in QSAT (due to 
    ! changes in T and heating effect of Q cond/evap) AND the change in Q)
    ! Note that condensation is -ve, evaporation is positive in this convention
    ZDQS(JL) = ZQSATNEW - ZQTMP(JL)

  ENDDO

  !-------------------------------------------------------
  !
  ! Condensation in clear air part of the grid box
  ! increasing cloud cover
  !
  !-------------------------------------------------------
  
  DO JL=KIDIA,KFDIA

    ZLCOND_ENV_L(JL) = 0.0_JPRB
    ZLCOND_ENV_I(JL) = 0.0_JPRB
    ZACOND(JL)       = 0.0_JPRB
    
    IF(ZDQS(JL) <= -RLMIN .AND. PA(JL,JK)<1.0_JPRB-ZEPSEC) THEN


      !---------------------------
      ! Supersaturation options
      !---------------------------      
      ZQE = (PQ(JL,JK)-PA(JL,JK)*ZQSMIX(JL))/&
          & MAX(ZEPSEC,1.0_JPRB-PA(JL,JK))  
      ZQE = MAX(0.0_JPRB,ZQE)
      
      IF(ZQE >= ZRHC(JL)*ZQSLIQK(JL) .AND. ZQE<ZQSLIQK(JL)) THEN
        ! note: not **2 on 1-a term if ZQE is used. 
        ZACOND(JL) = -(1.0_JPRB-PA(JL,JK))*ZDQS(JL)/&
         &MAX(2.0_JPRB*(ZQSLIQK(JL)-ZQE),ZEPSEC)

        ZACOND(JL) = MIN(ZACOND(JL),1.0_JPRB-PA(JL,JK))

        ZLCOND_ENV(JL) = -ZDQS(JL)*0.5_JPRB*ZACOND(JL)

        ! new limiter formulation
        ZZDL = 2.0_JPRB*(ZQSLIQK(JL)-ZQE)/MAX(ZEPSEC,1.0_JPRB-PA(JL,JK))
        IF (ZDQS(JL) < -ZZDL) THEN
          ZLCONDLIM = (PA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSLIQK(JL)+PQ(JL,JK)
          ZLCOND_ENV(JL) = MIN(ZLCOND_ENV(JL),ZLCONDLIM)
        ENDIF
        ZLCOND_ENV(JL) = MAX(ZLCOND_ENV(JL),0.0_JPRB)

        IF(ZLCOND_ENV(JL) < RLMIN .OR. (1.0_JPRB-PA(JL,JK))<ZEPSEC ) THEN
          ZLCOND_ENV(JL) = 0.0_JPRB
          ZACOND(JL)     = 0.0_JPRB
        ENDIF
        IF(ZLCOND_ENV(JL) == 0.0_JPRB) ZACOND(JL)=0.0_JPRB

        !------------------------------------------------------------------------
        ! All increase goes into liquid unless so cold cloud homogeneously freezes
        ! Include new liquid formation in first guess value, otherwise liquid 
        ! remains at cold temperatures until next timestep.
        !------------------------------------------------------------------------
        IF (PT(JL,JK)>RTHOMO) THEN
          ZLCOND_ENV_L(JL) = ZLCOND_ENV(JL)
        ELSE
          ZLCOND_ENV_I(JL) = ZLCOND_ENV(JL)
        ENDIF

      ENDIF
    ENDIF

  ENDDO

  !---------------------------------------------------------------------------
  !
  !  Condensation in cloudy part of the gridbox if ds < 0
  !
  !---------------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    
    ! Initialise output arrays
    ZLCOND_CLD_L(JL) = 0.0_JPRB
    ZLCOND_CLD_I(JL) = 0.0_JPRB

    IF(PA(JL,JK) > ZEPSEC .AND. ZDQS(JL) <= -RLMIN) THEN

      ZLCOND_CLD(JL) = MAX(-ZDQS(JL),0.0_JPRB)
   
      ! This is the maximum condensation that could occur
      ZCDMAX = MAX((ZQ(JL)-ZQSMIX_UPD(JL))/ZCORQSMIX_UPD(JL),0.0_JPRB)

      ZLCOND_CLD(JL) = MIN(ZLCOND_CLD(JL),ZCDMAX)

      ZLCOND_CLD(JL) = PA(JL,JK)*ZLCOND_CLD(JL)
      
      IF (ZLCOND_CLD(JL) < RLMIN) ZLCOND_CLD(JL)=0.0_JPRB

      !-------------------------------------------------------------------------
      ! All increase goes into liquid unless so cold cloud homogeneously freezes
      ! Include new liquid formation in first guess value, otherwise liquid 
      ! remains at cold temperatures until next timestep.
      !-------------------------------------------------------------------------
      IF (PT(JL,JK)>RTHOMO) THEN
        ZLCOND_CLD_L(JL) = ZLCOND_CLD(JL)
      ELSE
        ZLCOND_CLD_I(JL) = ZLCOND_CLD(JL)
      ENDIF

    ENDIF

  ENDDO

  !---------------------------------------------------------------------------
  !
  ! Evaporation of cloudy part of the gridbox if ds > 0
  !
  !---------------------------------------------------------------------------

  DO JL=KIDIA,KFDIA

    IF (ZDQS(JL) > 0.0_JPRB .AND. ZLI(JL) > 0.0_JPRB) THEN

      ZLEVAP_CLD(JL)   = MIN(PA(JL,JK)*ZDQS(JL),ZLI(JL))
      ZLEVAP_CLD_L(JL) = ZLIQFRAC(JL)*ZLEVAP_CLD(JL)
      ZLEVAP_CLD_I(JL) = ZICEFRAC(JL)*ZLEVAP_CLD(JL)

      ! Assumes a delta function of cloud condensate, no change to cloud cover
      ! If all condensate evaporated, then remove all cloud cover
      ZAEVAP(JL) = 0.0_JPRB

    ELSE
    
      ZLEVAP_CLD_L(JL) = 0.0_JPRB
      ZLEVAP_CLD_I(JL) = 0.0_JPRB
      ZAEVAP(JL)       = 0.0_JPRB
    
    ENDIF

  ENDDO



  !======================================================================
  !
  !
  !  3.  UPDATE VARIABLES FOR SUPERSATURATION CHECK
  !
  !
  !======================================================================

  !------------------------------------------------------------------
  ! Add tendencies for supersaturation check
  !------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
  
    ! Update temperature
    ZT_UPD(JL) = ZT(JL) &
              & + RALVDCP*(ZLCOND_CLD_L(JL)+ZLCOND_ENV_L(JL)-ZLEVAP_CLD_L(JL)) &
              & + RALSDCP*(ZLCOND_CLD_I(JL)+ZLCOND_ENV_I(JL)-ZLEVAP_CLD_I(JL))

    ! Update humidity
    ZQ_UPD(JL) = ZQ(JL) &
              & + (ZLEVAP_CLD_L(JL)-ZLCOND_CLD_L(JL)-ZLCOND_ENV_L(JL)) &
              & + (ZLEVAP_CLD_I(JL)-ZLCOND_CLD_I(JL)-ZLCOND_ENV_I(JL))

    ! Update cloud cover (ensure between bounds of 0 and 1)
    ZANEW = ZA(JL)+ZACOND(JL)-ZAEVAP(JL)
    ZANEW = MIN(ZANEW,1.0_JPRB)
    IF (ZANEW<RAMIN) ZANEW=0.0_JPRB
    ZA_UPD(JL) = ZANEW
  
  ENDDO
  
  !======================================================================
  !
  !
  !  4.  SUPERSATURATION CHECK
  !
  !
  !======================================================================

  ! Saturation adjustment step if there is any supersaturation above the defined limit
  IFTLIQICE = 2  ! Distribute liquid and ice, all liquid warmer than homog freezing temperature
  CALL CLOUD_SUPERSATCHECK(YDECLDP, KIDIA, KFDIA, KLON, KLEV, IFTLIQICE, &
                         & ZT_UPD, ZQ_UPD, ZA_UPD, PAP(:,JK), PAPH(:,KLEV+1), &
                         & ZT_ADJ, ZQ_ADJ, ZA_ADJ, ZL_ADJ, ZI_ADJ)
                             
  DO JL=KIDIA,KFDIA
    
    !-------------------------------------------------------------------
    ! Here the supersaturation is turned into liquid water
    ! However, if the temperature is below the threshold for homogeneous
    ! freezing then the supersaturation is turned instantly to ice.
    !--------------------------------------------------------------------
    IF (ZT_UPD(JL) > RTHOMO) THEN
      ! Turn supersaturation into liquid water        
      ZSUPSATL(JL) = ZL_ADJ(JL) + ZI_ADJ(JL)
      ZSUPSATI(JL) = 0.0_JPRB
    ELSE
      ! Turn supersaturation into ice water        
      ZSUPSATL(JL) = 0.0_JPRB
      ZSUPSATI(JL) = ZL_ADJ(JL) + ZI_ADJ(JL)
    ENDIF

    ! Calculate increase in cloud amount
    ZSUPSATA(JL) = ZA_ADJ(JL)

  ENDDO ! on JL


  !------------------------------------------------------------------
  ! Calculate total tendencies
  !------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
  
    ! Temperature tendency
    PTENDENCY_LOC_T(JL,JK) = &
        &   RALVDCP*(ZLCOND_CLD_L(JL)+ZLCOND_ENV_L(JL)+ZSUPSATL(JL)-ZLEVAP_CLD_L(JL))*ZQTMST &
        & + RALSDCP*(ZLCOND_CLD_I(JL)+ZLCOND_ENV_I(JL)+ZSUPSATI(JL)-ZLEVAP_CLD_I(JL))*ZQTMST
  
    ! Humidity tendency
    PTENDENCY_LOC_Q(JL,JK) = &
        &  (ZLEVAP_CLD_L(JL)-ZLCOND_CLD_L(JL)-ZLCOND_ENV_L(JL)-ZSUPSATL(JL) &
        & + ZLEVAP_CLD_I(JL)-ZLCOND_CLD_I(JL)-ZLCOND_ENV_I(JL)-ZSUPSATI(JL))*ZQTMST

    ! Cloud cover (ensure between bounds of 0 and 1)
    ZANEW = ZA(JL)+ZACOND(JL)-ZAEVAP(JL)+ZSUPSATA(JL)
    ZANEW = MIN(ZANEW,1.0_JPRB)
    IF (ZANEW<RAMIN) ZANEW=0.0_JPRB
    ZDA(JL)=ZANEW-ZA(JL)
    PTENDENCY_LOC_A(JL,JK) = ZDA(JL)*ZQTMST

    ! Cloud liquid and ice phase
    PTENDENCY_LOC_L(JL,JK) = (ZLCOND_CLD_L(JL)+ZLCOND_ENV_L(JL)+ZSUPSATL(JL) &
        & - ZLEVAP_CLD_L(JL))*ZQTMST
    PTENDENCY_LOC_I(JL,JK) = (ZLCOND_CLD_I(JL)+ZLCOND_ENV_I(JL)+ZSUPSATI(JL) &
        & - ZLEVAP_CLD_I(JL))*ZQTMST

    ! Tidy up any small negative values due to numerical truncation
    !Cloud liquid
    ZL_ADJ(JL) = ZL(JL)+PTENDENCY_LOC_L(JL,JK)*PTSPHY
    IF (ZL_ADJ(JL) < -RLMIN) WRITE(NULOUT,*) 'CLOUD_SATADJ: WARNING Negative cloud liquid! ',JK,ZL_ADJ(JL)
    IF (ZL_ADJ(JL) < 0.0_JPRB) THEN
      PTENDENCY_LOC_L(JL,JK) = PTENDENCY_LOC_L(JL,JK) - ZL_ADJ(JL)
      PTENDENCY_LOC_Q(JL,JK) = PTENDENCY_LOC_Q(JL,JK) + ZL_ADJ(JL)
      PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK) - RALVDCP*ZL_ADJ(JL)
    ENDIF
    ! Cloud ice
    ZI_ADJ(JL) = ZI(JL)+PTENDENCY_LOC_I(JL,JK)*PTSPHY
    IF (ZI_ADJ(JL) < -RLMIN) WRITE(NULOUT,*) 'CLOUD_SATADJ: WARNING Negative cloud ice! ',JK,ZI_ADJ(JL)
    IF (ZI_ADJ(JL) < 0.0_JPRB) THEN
      PTENDENCY_LOC_I(JL,JK) = PTENDENCY_LOC_I(JL,JK) - ZI_ADJ(JL)
      PTENDENCY_LOC_Q(JL,JK) = PTENDENCY_LOC_Q(JL,JK) + ZI_ADJ(JL)
      PTENDENCY_LOC_T(JL,JK) = PTENDENCY_LOC_T(JL,JK) - RALSDCP*ZI_ADJ(JL)
    ENDIF
 
  ENDDO

  !-----------------------------------------------------------------
  !
  !                   Diagnostics for cloud budget
  !           (currently hardwired, will be generalised)
  !
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Vertical integral of cloud process terms in one 3D field 
  !-----------------------------------------------------------------

  ! set pointer in extra variables array
  IF (LBUD23) THEN
    IS = 26  ! Take account of LBUD23 array diagnostics if turned on
  ELSE
    IS = 0
  ENDIF
  
  IF (LCLDBUD_VERTINT) THEN
 
   DO JL=KIDIA,KFDIA
  
    ! Cloud fraction budget terms
    IK = 0
    PEXTRA(JL,IK+1,1)  = PEXTRA(JL,IK+1,1)  + ZACOND(JL)*ZQTMST*ZDZ(JL)   ! + Condensation of new cloud
    PEXTRA(JL,IK+2,1)  = PEXTRA(JL,IK+2,1)  + ZAEVAP(JL)*ZQTMST*ZDZ(JL)   ! - Evaporation of cloud
    PEXTRA(JL,IK+3,1)  = PEXTRA(JL,IK+3,1)  + ZSUPSATA(JL)*ZQTMST*ZDZ(JL) ! + Supersat clipping after cloud_satadj
       
    ! Cloud condensate budget terms
    IK = 14
    PEXTRA(JL,IK+1,1) = PEXTRA(JL,IK+1,1) + ZLCOND_ENV_L(JL)*ZQTMST*ZDZ(JL) ! + Condensation of new cloud
    PEXTRA(JL,IK+2,1) = PEXTRA(JL,IK+2,1) + ZLCOND_CLD_L(JL)*ZQTMST*ZDZ(JL) ! + Condensation of existing cloud
    PEXTRA(JL,IK+3,1) = PEXTRA(JL,IK+3,1) - ZLEVAP_CLD_L(JL)*ZQTMST*ZDZ(JL)  ! - Evaporation of existing cloud
    PEXTRA(JL,IK+4,1) = PEXTRA(JL,IK+4,1) + ZSUPSATL(JL)*ZQTMST*ZDZ(JL) ! + Supersat clipping after cloud_satadj
    
    IK = 39
    PEXTRA(JL,IK+1,1) = PEXTRA(JL,IK+1,1) + ZLCOND_ENV_I(JL)*ZQTMST*ZDZ(JL) ! + Condensation of new cloud
    PEXTRA(JL,IK+2,1) = PEXTRA(JL,IK+2,1) + ZLCOND_CLD_I(JL)*ZQTMST*ZDZ(JL) ! + Condensation of existing cloud
    PEXTRA(JL,IK+3,1) = PEXTRA(JL,IK+3,1) - ZLEVAP_CLD_I(JL)*ZQTMST*ZDZ(JL)  ! - Evaporation of existing cloud
    PEXTRA(JL,IK+4,1) = PEXTRA(JL,IK+4,1) + ZSUPSATI(JL)*ZQTMST*ZDZ(JL) ! + Supersat clipping after cloud_satadj

   ENDDO
   IS = IS + 1

  ENDIF ! on LCLDBUD_VERTINT

  !-----------------------------------------------------------------
  ! Cloud fraction budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDC) THEN
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+1)  = PEXTRA(JL,JK,IS+1) + ZACOND(JL)*ZQTMST   ! Condensation of new cloud
      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZAEVAP(JL)*ZQTMST   ! Evaporation of cloud
      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZSUPSATA(JL)*ZQTMST ! Supersat clipping
      !PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZDTDP*ZWTOT
      !PEXTRA(JL,JK,IS+14) = PEXTRA(JL,JK,IS+14) + ZDTDIAB*ZQTMST
      !PEXTRA(JL,JK,IS+15) = PEXTRA(JL,JK,IS+15) + PTENDENCY_VDF_T(JL,JK)
      !PEXTRA(JL,JK,IS+16) = PEXTRA(JL,JK,IS+16) + PLUDELI(JL,JK,4)*ZGDP(JL)
      !PEXTRA(JL,JK,IS+17) = PEXTRA(JL,JK,IS+17) + PTENDENCY_VDF_Q(JL,JK)
      !PEXTRA(JL,JK,IS+18) = PEXTRA(JL,JK,IS+18) + PLUDELI(JL,JK,3)*ZGDP(JL)
    ENDDO
    IS = IS + 12
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud liquid condensate budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDL) THEN
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+1) = PEXTRA(JL,JK,IS+1) + ZLCOND_ENV_L(JL)*ZQTMST ! + Condensation of new cloud
      PEXTRA(JL,JK,IS+2) = PEXTRA(JL,JK,IS+2) + ZLCOND_CLD_L(JL)*ZQTMST ! + Condensation of existing cloud
      PEXTRA(JL,JK,IS+3) = PEXTRA(JL,JK,IS+3) - ZLEVAP_CLD_L(JL)*ZQTMST ! - Evaporation of existing cloud
      PEXTRA(JL,JK,IS+4) = PEXTRA(JL,JK,IS+4) + ZSUPSATL(JL)*ZQTMST     ! + Supersat clipping so far this timestep
   ENDDO
    IS = IS + 22
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud ice condensate budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDI) THEN
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+1)  = PEXTRA(JL,JK,IS+1) + ZLCOND_ENV_I(JL)*ZQTMST ! + Condensation of new cloud
      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZLCOND_CLD_I(JL)*ZQTMST ! + Condensation of existing cloud
      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) - ZLEVAP_CLD_I(JL)*ZQTMST ! - Evaporation of existing cloud
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + ZSUPSATI(JL)*ZQTMST     ! + Supersat clipping so far this timestep 
    ENDDO
    IS = IS + 18
  ENDIF
 
ENDDO ! on vertical level JK

!===============================================================================
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ',1,ZHOOK_HANDLE)
!===============================================================================
CONTAINS
FUNCTION QSATMIXADJ(PT,PQ,PP)

!------------------------------------------------------------------------------- 
! Description:
! 
! Performs saturation adjustment with respect to mixed-phase function
!
! Inputs: PT - temperature (K)
!         PQ - specific humidity (kg/kg)
!         PP - pressure (Pa)
! 
! Returns the amount of condensed water to bring the air to water saturation
! 
! Author:
!   R. Forbes Aug 2018
!
!------------------------------------------------------------------------------- 

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RETV
USE YOETHF    ,ONLY : R4LES, R5LES, RALVDCP

IMPLICIT NONE

REAL(KIND=JPRB),INTENT(IN) :: PT
REAL(KIND=JPRB),INTENT(IN) :: PQ
REAL(KIND=JPRB),INTENT(IN) :: PP

REAL(KIND=JPRB) :: QSATMIXADJ
REAL(KIND=JPRB) :: ZT, ZQ, ZP_R, ZFOEET, ZCOR, ZQSAT, ZFACW, ZDQSDT, ZCORQS, ZCOND

INTEGER(KIND=JPIM) :: JI, I_NITER_SATADJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------- 
!IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ:QSATMIXADJ',0,ZHOOK_HANDLE)

! Number of saturation adjustment iterations 
I_NITER_SATADJ = 2

ZT  = PT
ZQ  = PQ
ZP_R = 1.0_JPRB/PP  
 
! Loop over number of iterations
DO JI = 1,I_NITER_SATADJ 
  
  ! Calculate saturation from updated temperature
  ZFOEET = MIN(FOEEWM(ZT)*ZP_R,0.5_JPRB)
  ZCOR   = 1.0_JPRB/(1.0_JPRB - RETV*ZFOEET)
  ZQSAT  = ZFOEET*ZCOR

  ! Saturation adjustment term to bring gridbox T,Q to saturation
  ZCORQS = 1.0_JPRB+ZQSAT*ZCOR*FOEDEM(ZT)
  ! Calculate amount of condensed water from first iteration
  ZCOND  = (ZQ-ZQSAT)/ZCORQS

  ! Calculate new values of temperature and humidity
  ZT = ZT + FOELDCPM(ZT)*ZCOND
  ZQ = ZQ - ZCOND
 
ENDDO

! Set output of function to change in humidity/condensed water
QSATMIXADJ = ZQ  

!IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ:QSATMIXADJ',1,ZHOOK_HANDLE)
END FUNCTION QSATMIXADJ
!===============================================================================
FUNCTION QSATWATADJ(PT,PQ,PP)

!------------------------------------------------------------------------------- 
! Description:
! 
! Performs saturation adjustment with respect to water
!
! Inputs: PT - temperature (K)
!         PQ - specific humidity (kg/kg)
!         PP - pressure (Pa)
! 
! Returns the amount of condensed water to bring the air to water saturation
! 
! Author:
!   R. Forbes Aug 2018
!
!------------------------------------------------------------------------------- 

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RETV
USE YOETHF    ,ONLY : R4LES, R5LES, RALVDCP

IMPLICIT NONE

REAL(KIND=JPRB),INTENT(IN) :: PT
REAL(KIND=JPRB),INTENT(IN) :: PQ
REAL(KIND=JPRB),INTENT(IN) :: PP

REAL(KIND=JPRB) :: QSATWATADJ
REAL(KIND=JPRB) :: ZT, ZQ, ZP_R, ZFOEET, ZCOR, ZQSAT, ZFACW, ZDQSDT, ZCORQS, ZCOND

INTEGER(KIND=JPIM) :: JI, I_NITER_SATADJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------- 
!IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ:QSATWATADJ',0,ZHOOK_HANDLE)

! Number of saturation adjustment iterations 
I_NITER_SATADJ = 2

ZT  = PT
ZQ  = PQ
ZP_R = 1.0_JPRB/PP  
 
! Loop over number of iterations
DO JI = 1,I_NITER_SATADJ 
  
  ! Calculate saturation from updated temperature
  ZFOEET = MIN(FOEELIQ(ZT)*ZP_R,0.5_JPRB)
  ZCOR   = 1.0_JPRB/(1.0_JPRB - RETV*ZFOEET)
  ZQSAT  = ZFOEET*ZCOR

  ! Saturation adjustment term to bring gridbox T,Q to saturation
  ZFAC   = R5LES/((ZT-R4LES)**2)
  ZDQSDT = ZFAC*ZCOR*ZQSAT
  ZCORQS = 1.0_JPRB + RALVDCP*ZDQSDT
  ! Calculate amount of condensed water from first iteration
  ZCOND  = (ZQ-ZQSAT)/ZCORQS

  ! Calculate new values of temperature and humidity
  ZT = ZT + RALVDCP*ZCOND
  ZQ = ZQ - ZCOND
 
ENDDO

! Set output of function to change in humidity/condensed water
QSATWATADJ = ZQ  

!IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ:QSATWATADJ',1,ZHOOK_HANDLE)
END FUNCTION QSATWATADJ

END SUBROUTINE CLOUD_SATADJ
