! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE AER_DRYDEPVELZH14(KSEASON_WE,KVEG,KVEGZH,PLAI, PWETD, PSIGMA, PRHO,&
      &    PZ0M,PCI,PUSTR,PT,PDZ,PDUST_REBOUND,PVDEP)


!**** *AER_DRYDEPVELZH14* -  ROUTINE FOR PARAMETRIZATION OF DRY DEPOSITION VELOCITY

!**   DESCRIPTION 
!     ----------
! Calculates aerosol dry deposition (and sedimentation).
! Based on the parameterisation of Zhang and He (2014)
!
!**   INTERFACE.
!     ----------
!          *AER_DRYDEPVELZH14* IS CALLED FROM *AER_DRYDEP*.

! INPUTS:
! -------
! KSEASON    : Season
! KVEG       : IFS land class  (index)
! KVEGZH     : Zhang01 land class  (index)
! PWETD      : PARTICLE WET DIAMETER   (m)
! PSIGMA     : Sigmag
! PRHOi      : Particle density (kg/m3)
! PLAIi      : Leaf Area Index
! PZ0M       : ROUGHNESS LENGTH         (m)
! PCI        : SEA-ICE MASK
! PT         : Temperature (K)
! PUST       : FRICTION VELOCITY        (m.s-1)
! PDZ        : DELTA Z                  (m)
! PDUST_REBOUND        : Fraction of super coarse dust that is subjected to rebound effect

! OUTPUTS:
! --------
! PVDEP      : DRY DEPOSITION VELOCITY  (m.s-1)

! LOCAL VARIABLES:
! ---------------
! ZSR_AV_3    : 3rd moment avg surface resistance
! ZAR         : Aerodynamic resistance
!


!     AUTHOR.
!     -------
!        SAMUEL REMY   *HYGEOS*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2019-06-28

! 
!     REFERENCES.
!     ----------
! Zhang and He, ACP, 2014
!
!----------------------------------------------------------------------
USE PARKIND1,            ONLY: JPRB, JPRD, JPIM
USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK
USE YOMLUN,              ONLY: NULOUT
USE YOMCST ,             ONLY : RG, RPI, RKBOL

IMPLICIT NONE


! .. Subroutine interface
INTEGER(KIND=JPIM), INTENT(IN) :: KSEASON_WE
INTEGER(KIND=JPIM), INTENT(IN) :: KVEG
INTEGER(KIND=JPIM), INTENT(IN) :: KVEGZH
REAL(KIND=JPRB), INTENT(IN) :: PLAI
REAL(KIND=JPRB), INTENT(IN) :: PWETD
REAL(KIND=JPRB), INTENT(IN) :: PSIGMA
REAL(KIND=JPRB), INTENT(IN) :: PRHO
REAL(KIND=JPRB), INTENT(IN) :: PZ0M
REAL(KIND=JPRB), INTENT(IN) :: PCI
REAL(KIND=JPRB), INTENT(IN) :: PUSTR
REAL(KIND=JPRB), INTENT(IN) :: PT
REAL(KIND=JPRB), INTENT(IN) :: PDZ
REAL(KIND=JPRB), INTENT(IN) :: PDUST_REBOUND
REAL(KIND=JPRB), INTENT(OUT) :: PVDEP
!
!    Local Variables
REAL(KIND=JPRB)    :: ZKARMN
REAL(KIND=JPRB)    :: ZAR, ZVDS, ZK1
REAL(KIND=JPRB)    :: ZSR_AV_3,ZVISC, ZMFPA, ZVBA,ZVGRAV_AV_3,ZSN_AV_3,ZREBOUND
REAL(KIND=JPRB)    :: ZCR(15,5)
REAL(KIND=JPRB), PARAMETER :: ZIMA=4.78E-26_JPRB
REAL(KIND=JPRB), PARAMETER :: ZA1(5) = (/3.4E-3_JPRB, 4.3E-3_JPRB, 4.8E-3_JPRB,5.4E-3_JPRB, 6.9E-3_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZB1(15) = (/2.6E-1_JPRB, 3.9E-1_JPRB, -1.6E-1_JPRB,1.6E-2_JPRB, -5.3E-2_JPRB,&
                                         &6.7E-2_JPRB,5.6E-2_JPRB,7.5E-2_JPRB,7.1E-2_JPRB, 9.9E-2_JPRB, &
                                         &-1.2E-1_JPRB, 1.6E-2_JPRB,5.6E-2_JPRB, -7.9E-2_JPRB,-6.E-2_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZB2(15) = (/-1.3_JPRB, -3.3_JPRB,1.5_JPRB, 3.4E-1_JPRB, 6.6E-1_JPRB, &
                                         &3.2E-2_JPRB,1.6E-1_JPRB,1.2E-1_JPRB,7.0E-3_JPRB,-1.3E-2_JPRB,&
                                         &1.2_JPRB,3.4E-1_JPRB,1.6E-1_JPRB,1.0_JPRB,1.0_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZB3(15) = (/3.0_JPRB, 8.8_JPRB,7.8E-1_JPRB, 4.5E-1_JPRB, 6.7E-1_JPRB, &
                                         &1.2E-1_JPRB,2.8E-1_JPRB,2.4E-1_JPRB,5.7E-2_JPRB,4.6E-2_JPRB,&
                                         &7.1E-1_JPRB,4.5E-1_JPRB,2.8E-1_JPRB,6.6E-1_JPRB,6.5E-1_JPRB/)

REAL(KIND=JPRB), PARAMETER :: ZD1(15) = (/-8.7E-1_JPRB, -7.3_JPRB, -9.8E-1_JPRB,-2.2_JPRB, -1.7_JPRB,&
                                         &-1.3_JPRB,-2.2_JPRB,-2.1_JPRB,-7.2E-1_JPRB, -9.8E-2_JPRB, &
                                         &-1.6_JPRB, -2.2_JPRB,-2.2_JPRB, -2.0_JPRB,-2.0_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZD2(15) = (/-5.5_JPRB, 4.6E1_JPRB,7.1E1_JPRB, 3.9E1_JPRB, 5.2E1_JPRB, &
                                         &1.3E1_JPRB,2.7E1_JPRB,2.4E1_JPRB,6.4_JPRB,2.1_JPRB,&
                                         &6.6E1_JPRB,3.9E1_JPRB,2.7E1_JPRB,6.3E1_JPRB,6.2E1_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZD3(15) = (/9.9E1_JPRB, 9.4E1_JPRB,-9.5_JPRB, -6.7_JPRB, -1.2E1_JPRB, &
                                         &5.3E-1_JPRB,-2.7_JPRB,-1.8_JPRB,1.4_JPRB,3.3_JPRB,&
                                         &-1.7E1_JPRB,-6.7_JPRB,-2.7_JPRB,-1.6E1_JPRB,-1.5E1_JPRB/)

REAL(KIND=JPRB), PARAMETER :: ZC1(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &4.8_JPRB, 1.8_JPRB,7.4E-1_JPRB, 5.1_JPRB,3.4_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZC2(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &-5.1_JPRB, -2.0E-1_JPRB,1.7_JPRB, -4.2_JPRB,-2.4_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZC3(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &1.8_JPRB, -5.3E-1_JPRB,1.4_JPRB, 9.9E-1_JPRB,3.4E-1_JPRB/)

REAL(KIND=JPRB), PARAMETER :: ZF1(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &7.7_JPRB, 6.2_JPRB,7.7E-1_JPRB, 1.1E1_JPRB,7.9_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZF2(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &-1.5E1_JPRB, -1.2E1_JPRB,-1.4E1_JPRB, 2.0E1_JPRB,-1.5E1_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZF3(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &7.8_JPRB, 6.1E-1_JPRB,7.4_JPRB, 1.1E1_JPRB,8.0_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZLAIMAX(15) = (/0._JPRB, 0._JPRB, 0._JPRB,0._JPRB, 0._JPRB,&
                                         &0._JPRB, 0._JPRB,0._JPRB,0._JPRB, 0._JPRB, &
                                         &5._JPRB, 5._JPRB,3._JPRB, 2._JPRB,4._JPRB/)

INTEGER(KIND=JPIM) :: JLAND
INTEGER(KIND=JPIM) :: JLAND2

REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "aer_vgrav.intfb.h"


IF (LHOOK) CALL DR_HOOK('AER_DRYDEPVELZH14',0,ZHOOK_HANDLE)

! IFS Land classed
!1)=100._JPRB   L ! Crops, Mixed Farming             2: agricultural land
!2)=100._JPRB   L ! Short Grass                      3: range land
!3)=250._JPRB   H ! Evergreen Needleleaf Trees       5: coniferous forest
!4)=250._JPRB   H ! Deciduous Needleleaf Trees       5: coniferous forest
!5)=175._JPRB   H ! Deciduous Broadleaf Trees        4: deciduous forest
!6)=240._JPRB   H ! Evergreen Broadleaf Trees        4: deciduous forest
!7)=100._JPRB   L ! Tall Grass                       3: range land
!8)=250._JPRB    ! Desert
!9)=80._JPRB    L ! Tundra                          11: rocky open areas with
!low shrubs
!10)=100._oPRB  L ! Irrigated Crops                  2: agricultural
!11)=150._JPRB  L ! Semidesert                       8: barren land / desert
!12)=0.0_JPRB      ! Ice Caps and Glaciers
!13)=240._JPRB  L ! Bogs and Marshes                 9: nonforested wetland
!14)=0.0_JPRB      ! Inland Water
!15)=0.0_JPRB      ! Ocean
!16)=225._JPRB  L ! Evergreen Shrubs                 3: range land
!17)=225._JPRB  L ! Deciduous Shrubs                 3: range land
!18)=250._JPRB  H ! Mixed Forest/woodland            6: mixed forest & wet lands
!19)=175._JPRB  H ! Interrupted Forest               6: mixed forest & wet lands
!20)=150._JPRB  L ! Water and Land Mixtures          9: nonforested wetland

ZCR(:,1)=(/8E-4_JPRB,2e-3_JPRB,8E-4_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,2)=(/8E-4_JPRB,2e-3_JPRB,8E-4_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,3)=(/8E-4_JPRB,2e-3_JPRB,2E-3_JPRB,4E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,4)=(/8E-4_JPRB,2e-3_JPRB,2E-3_JPRB,4E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,5)=(/8E-4_JPRB,2e-3_JPRB,8E-4_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)


! map soil types
SELECT CASE(KVEG)
  CASE(12, 3, 6, 4, 5, 9, 8, 18, 19, 11)
    JLAND2=2
  CASE(16, 17, 7)
    JLAND2=3
  CASE(2, 1, 10, 13, 20)
    JLAND2=4
  CASE(15, 14)
    JLAND2=5
  CASE DEFAULT
    WRITE(NULOUT,*) 'AER_DRYDEP: unknown vegetation type KVEG=', KVEG
    CALL ABOR1('AER_DRYDEP: unknown vegetation type when setting JLAND2')
END SELECT

! map soil types
SELECT CASE(KVEG)
  CASE(14, 15) ! Ocean and inland lake
    JLAND=1
  CASE(8, 9, 12, 11) ! Ice cap Tundra desert semi-desert
    JLAND=2
  CASE(3) ! Evergreen needleleaf tree
    JLAND=3
  CASE(6, 18, 19) ! Evergreen broadleaf trees & Mixed woods & Interrupted forest
    JLAND=4
  CASE(4) ! Deciduous needleleaf trees
    JLAND=11
  CASE(5) ! Deciduous broadleaf trees
    JLAND=12
  CASE(1) ! Crops
    JLAND=15
  CASE(7) ! Long grass
    JLAND=14
  CASE(2, 10, 16) ! Short grass & Irrigated crops & evergeen shrubs
    JLAND=8
  CASE(17) ! Deciduous shrubs
    JLAND=11
  CASE(13, 20) ! Swamp & Water/Land mixtures
    JLAND=9
  CASE DEFAULT
    WRITE(NULOUT,*) 'AER_DRYDEP: unknown vegetation type KVEG=', KVEG
    CALL ABOR1('AER_DRYDEP: unknown vegetation type when setting JLAND')
END SELECT

PVDEP=0.0_JPRB
ZKARMN=0.4_JPRB
! .. Calculate aerodynamic resistance
ZAR=LOG(PDZ/PZ0M)/(ZKARMN*PUSTR)
IF (PWETD < 2.5E-6_JPRB) THEN
  ZVDS=ZA1(JLAND2)*PUSTR
ELSEIF (PWETD >= 2.5E-6 .AND. PWETD < 10.E-6_JPRB) THEN
  IF (JLAND > 10) THEN
    ZK1=ZC1(JLAND)*PUSTR + ZC2(JLAND)*PUSTR**2 + ZC3(JLAND)*PUSTR**3
    IF (JPRB/=JPRD) THEN
      ZK1=MIN(ZK1,20._JPRB)
      ZK1=MAX(ZK1,-20._JPRB)
    ENDIF
    ZVDS=(ZB1(JLAND)*PUSTR + ZB2(JLAND)*PUSTR**2 + ZB3(JLAND)*PUSTR**3)*EXP(ZK1*(MIN(PLAI/ZLAIMAX(JLAND), 1._JPRB) - 1._JPRB))
  ELSE
    ZVDS=ZB1(JLAND)*PUSTR + ZB2(JLAND)*PUSTR**2 + ZB3(JLAND)*PUSTR**3
  ENDIF
ELSE
  IF (JLAND > 10) THEN
    ZK1=ZF1(JLAND)*PUSTR + ZF2(JLAND)*PUSTR**2 + ZF3(JLAND)*PUSTR**3
    IF (JPRB/=JPRD) THEN
      ZK1=MIN(ZK1,20._JPRB)
      ZK1=MAX(ZK1,-20._JPRB)
    ENDIF
    ZVDS=(ZD1(JLAND)*PUSTR + ZD2(JLAND)*PUSTR**2 + ZD3(JLAND)*PUSTR**3)*EXP(ZK1*(MIN(PLAI/ZLAIMAX(JLAND), 1._JPRB) - 1._JPRB))
  ELSE
    ZVDS=ZD1(JLAND)*PUSTR + ZD2(JLAND)*PUSTR**2 + ZD3(JLAND)*PUSTR**3
  ENDIF
ENDIF
!       Calculate surface resistance
IF (ZVDS > 1.E-6_JPRB) THEN
  ZSR_AV_3=1.0_JPRB/(ZVDS)
ELSE
  ZSR_AV_3=100000._JPRB
ENDIF

! dynamic viscosity of air (kg m^-1 s^-1)
ZVISC=1.83E-5_JPRB*(416.16_JPRB/(PT+120.0_JPRB))*(SQRT(PT/296.16_JPRB)**3.0_JPRB)
! mean free speed of air molecules
ZVBA  = SQRT(8.0_JPRB*RKBOL*PT)/(RPI*ZIMA)
 ! mean free path of air molecules
ZMFPA   = 2.0_JPRB*ZVISC/(PRHO*ZVBA)

!       Calculate 3rd moment avg. grav. settling velocities
CALL AER_VGRAV(3,PWETD,PSIGMA,ZVISC,ZMFPA,PRHO,ZVGRAV_AV_3)


IF (KVEGZH==8 .OR. KVEGZH==9 .OR. KVEGZH==12 .OR. KVEGZH==13 .OR.KVEGZH==14 ) THEN
!        Calculate stokes number for smooth surfaces
  ZSN_AV_3=ZVGRAV_AV_3*PUSTR*PUSTR/ZVISC
ELSE
!        Calculate stokes number for vegetated surfaces
  ZSN_AV_3=ZVGRAV_AV_3*PUSTR/(RG*ZCR(KVEGZH,KSEASON_WE))
ENDIF

ZREBOUND=1._JPRB
IF (PRHO > 2300._JPRB .AND. PUSTR > 0.335_JPRB .AND. PWETD >=2.5E-6_JPRB .AND. JLAND > 1) THEN
  ZREBOUND=EXP(-1_JPRB*PDUST_REBOUND*(ZSN_AV_3)**0.5_JPRB)
ENDIF

ZSR_AV_3=ZSR_AV_3/ZREBOUND
!       Calculate deposition velocity
PVDEP=1.0_JPRB/(ZAR+ZSR_AV_3)


!       Calculate deposition velocity
PVDEP=1.0_JPRB/(ZAR+ZSR_AV_3)

IF (LHOOK) CALL DR_HOOK('AER_DRYDEPVELZH14',1,ZHOOK_HANDLE)
END SUBROUTINE AER_DRYDEPVELZH14
