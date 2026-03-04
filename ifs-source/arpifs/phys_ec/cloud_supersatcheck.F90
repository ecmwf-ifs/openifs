! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLOUD_SUPERSATCHECK(YDECLDP, KIDIA, KFDIA, KLON, KLEV, KFTLIQICE, &
                             &  PT, PQ, PA, PP, PSURF, & 
                             &  PT_ADJ, PQ_ADJ, PA_ADJ, PL_ADJ, PI_ADJ)

!------------------------------------------------------------------------------- 
! Description:
! 
! Checks for supersaturation with respect to the defined gridbox mean 
! saturation limit and adjusts back to this saturation limit
!
! Inputs: PT - temperature (K)
!         PQ - specific humidity (kg/kg)
!         PA - cloud fraction (0-1)
!         PP - pressure (Pa)
! 
! Returns the changes to temperature, humidity, cloud fraction, 
!                        liquid and ice condensate
! 
! Author:
!   R. Forbes Jan 2020
!
!------------------------------------------------------------------------------- 

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOECLDP  , ONLY : TECLDP
USE YOMCST   , ONLY : RETV, RTT, RLVTT, RLSTT
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2

IMPLICIT NONE

! Input/output arguments
TYPE(TECLDP)      ,INTENT(IN) :: YDECLDP
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA        ! Start array location
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA        ! End array location
INTEGER(KIND=JPIM),INTENT(IN) :: KLON         ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV         ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN) :: KFTLIQICE    ! Controls T-dependence for liquid/ice production
REAL(KIND=JPRB),   INTENT(IN) :: PT(KLON)     ! Input temperature (K)
REAL(KIND=JPRB),   INTENT(IN) :: PQ(KLON)     ! Input humidity (K)
REAL(KIND=JPRB),   INTENT(IN) :: PA(KLON)     ! Input cloud fraction (0-1)
REAL(KIND=JPRB),   INTENT(IN) :: PP(KLON)     ! Input pressure (Pa)
REAL(KIND=JPRB),   INTENT(IN) :: PSURF(KLON)  ! Input surface pressure (Pa)
REAL(KIND=JPRB),   INTENT(OUT):: PT_ADJ(KLON) ! Output temperature change (K)
REAL(KIND=JPRB),   INTENT(OUT):: PQ_ADJ(KLON) ! Output humidity change (kg/kg)
REAL(KIND=JPRB),   INTENT(OUT):: PA_ADJ(KLON) ! Output cloud fraction change (0-1)
REAL(KIND=JPRB),   INTENT(OUT):: PL_ADJ(KLON) ! Output cloud water change (kg/kg)
REAL(KIND=JPRB),   INTENT(OUT):: PI_ADJ(KLON) ! Output cloud ice change (kg/kg)

! Local variables
REAL(KIND=JPRB) :: ZT        ! Gridbox mean temperature
REAL(KIND=JPRB) :: ZQ        ! Gridbox mean humidity
REAL(KIND=JPRB) :: ZA        ! Cloud fraction
REAL(KIND=JPRB) :: ZP_R      ! Reciprocal of atmospheric pressure
REAL(KIND=JPRB) :: ZQSLIQ    ! Saturation wrt water
REAL(KIND=JPRB) :: ZQSICE    ! Saturation wrt water (T>0), ice (T<0)
REAL(KIND=JPRB) :: ZQSLIQK   ! Saturation wrt to water (T>-40), Koop (T<-40)
REAL(KIND=JPRB) :: ZQSMIX    ! Saturation wrt water (T>0), mix ice/water, ice (T<RTICE)
REAL(KIND=JPRB) :: ZFOEELIQ  ! Factor for QSLIQ
REAL(KIND=JPRB) :: ZFOEEICE  ! Factor for QSICE
REAL(KIND=JPRB) :: ZFOKOOP   ! Factor for QSLIQK
REAL(KIND=JPRB) :: ZFOEALFA  ! Mixed phase function from water (=1) to ice (=0)
REAL(KIND=JPRB) :: ZFOEEMIX  ! Factor for QSMIX
REAL(KIND=JPRB) :: ZCOR      ! Factor for QSMIX
REAL(KIND=JPRB) :: ZCORQSMIX ! Condensation correction factor for QSMIX
!REAL(KIND=JPRB) :: ZCORQSLIQ ! Condensation correction factor for QSLIQ
REAL(KIND=JPRB) :: ZQSLIM    ! Saturation limit for partially cloudy gridbox 
REAL(KIND=JPRB) :: ZA_UPD    ! Updated cloud fraction if condensation occurs
REAL(KIND=JPRB) :: ZDA       ! Change in cloud fraction
REAL(KIND=JPRB) :: ZCOND     ! Condensation amount (kg/kg)
REAL(KIND=JPRB) :: ZFACI     ! Factor for ice cloud
!REAL(KIND=JPRB) :: ZFACW     ! Factor for water cloud
REAL(KIND=JPRB) :: ZSIGK     ! Sigma level
REAL(KIND=JPRB) :: ZRHC      ! Critical relative humidity
REAL(KIND=JPRB) :: ZEPSEC

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttre.func.h"

!------------------------------------------------------------------------------- 
IF (LHOOK) CALL DR_HOOK('CLOUD_SUPERSATCHECK',0,ZHOOK_HANDLE)

ASSOCIATE(RSSICEFACTOR=>YDECLDP%RSSICEFACTOR, RTHOMO=>YDECLDP%RTHOMO, &
        & RAMID=>YDECLDP%RAMID)

ZEPSEC = 1.E-14_JPRB

DO JL=KIDIA,KFDIA

  ZA_UPD = 0.0_JPRB
  ZDA    = 0.0_JPRB
 
  ZT   = PT(JL)
  ZQ   = PQ(JL)
  ZA   = PA(JL)
  ZP_R = 1.0_JPRB/PP(JL) 

  ! Safety checks
  ZT = MAX(ZT,160.0_JPRB)
  ZA = MAX(MIN(ZA,1.0_JPRB),0.0_JPRB)

  !---------------------------
  ! Critical relative humidity
  !---------------------------
  ZRHC     = RAMID
  ZSIGK    = PP(JL)/PSURF(JL)
  ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
  IF(ZSIGK > 0.8_JPRB) THEN
    ZRHC = RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
  ENDIF

  !-------------------------------------------------
  ! Calculate various saturation-related quantities
  !-------------------------------------------------
  
  ! Saturation wrt water (all T)
  ZFOEELIQ = MIN(FOEELIQ(ZT)*ZP_R,0.5_JPRB)
  ZQSLIQ   = ZFOEELIQ/(1.0_JPRB-RETV*ZFOEELIQ)

  ! Saturation wrt ice (T<0), water (T>0)
  ZFOEEICE = MIN(FOEEW(ZT)*ZP_R,0.5_JPRB)
  ZQSICE   = ZFOEEICE/(1.0_JPRB-RETV*ZFOEEICE)

  ! Saturation wrt Koop curve (T<0) liquid water saturation limited to Koop curve
  ! FOEDELTA = 1 for T>0, =0 for T<0
  ZFOKOOP = MIN(ZFOEEICE*(RKOOP1-RKOOP2*ZT),ZFOEELIQ)
  ZQSLIQK = ZFOKOOP/(1.0_JPRB-RETV*ZFOKOOP)

  ! Liquid-phase saturation adjustment factor to bring gridbox T,Q to saturation
  !ZFACW       = R5LES/((ZT-R4LES)**2)
  !ZCOR        = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEELIQ)
  !ZCORQSLIQ   = 1.0_JPRB+RALVDCP*ZFACW*ZCOR*ZQSLIQ

  ! Saturation wrt water (T>0), mixed (RTT<T<0), ice (T<RTT)
  ZFOEALFA = FOEALFA(ZT)
  ZFOEEMIX = MIN(FOEEWM(ZT)*ZP_R,0.5_JPRB)
  ZCOR     = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEMIX)
  ZQSMIX   = ZFOEEMIX*ZCOR

  ! Mixed-phase saturation adjustment factor to bring gridbox T,Q to saturation
  ZCORQSMIX = 1.0_JPRB+ZQSMIX*ZCOR*FOEDEM(ZT)

  ! Reduce maximum supersaturation to be half way between QSLIQ and QSMIX
  ! to reduce humidity in the upper
  ! Need to revisit this 
  ZQSLIQK = 0.5_JPRB*(ZQSLIQK+ZQSMIX)

  ! Calculate the maximum grid-box mean saturated state taking into account
  ! the saturation assumption in-cloud and ice supersaturated Koop limit
  ! in the clear air. Note ZQSLIQK = ZQSMIX for T warmer than 0degC
  ZQSLIM = ZA*ZQSMIX + (1.0_JPRB-ZA)*ZQSLIQK

  !-------------------------------------------
  ! If supersaturated, remove supersaturation
  !-------------------------------------------
  
  ! Determine if gridbox mean humidity is supersaturated
  ZCOND = MAX(ZQ-ZQSLIM,0.0_JPRB)

  IF (ZCOND > ZEPSEC) THEN 

    ! Recalculate cloud cover and new gridbox mean saturation limit,
    !  if T>0, cloud fraction must be 1 if grid-mean is supersaturated 
    !  if T<0, cloud fraction is reduced with timescale RKOOPTAU, i.e.
    !          the cloudy part is reduced to qsmix and the 
    !          clear air environment is reduced to qsliqk
    ZFACI = ZFOEALFA + (1.0_JPRB-ZFOEALFA)*RSSICEFACTOR
    ZDA   = (1.0_JPRB-ZA)*ZFACI
    
    ! Increase cloud amount using RKOOPTAU timescale
    ZA_UPD = ZA + ZDA

    ! Recalculate gridbox mean saturation limit
    ZQSLIM = ZA_UPD*ZQSMIX + (1.0_JPRB-ZA_UPD)*ZQSLIQK

    ! Calculate amount of condensation that can occur to leave the gridbox at QSLIM
    ZCOND = (ZQ-ZQSLIM)/ZCORQSMIX

  ENDIF

  !-----------------------------------------------------------
  ! Calculate new values of temperature, humidity, condensate
  !-----------------------------------------------------------
  
  PT_ADJ(JL) = FOELDCPM(ZT)*ZCOND
  PQ_ADJ(JL) = -ZCOND
  PA_ADJ(JL) = ZDA
  IF (KFTLIQICE == 1) THEN
    ! Partition liquid and ice according to alfa function 
    PL_ADJ(JL) = ZCOND*ZFOEALFA
    PI_ADJ(JL) = ZCOND*(1.0_JPRB-ZFOEALFA) 
  ELSEIF (KFTLIQICE == 2) THEN
    ! All liquid up to homogeneous freezing threshold, ice at colder temperatures
    IF (PT(JL)+PT_ADJ(JL) >= RTHOMO) THEN
      PL_ADJ(JL) = ZCOND
      PI_ADJ(JL) = 0.0_JPRB
    ELSE
      PL_ADJ(JL) = 0.0_JPRB
      PI_ADJ(JL) = ZCOND
    ENDIF 
  ENDIF

ENDDO ! on JL

!===============================================================================
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUD_SUPERSATCHECK',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUD_SUPERSATCHECK
