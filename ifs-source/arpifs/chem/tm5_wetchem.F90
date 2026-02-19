! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_WETCHEM(YGFL,KIDIA,KFDIA,KLON,PDT,PTEMP,PAP,PRS,PLP,PY0,PHPLUS,PY)
!**********************************************************************
!     
!TM5_wetchem - aqueous phase chemistry of sulfur  (and other)
!programmed by Ad Jeuken (KNMI) and Frank Dentener (IMAU)
!adapted for TM5 by Maarten Krol (IMAU) 1-2002
!
!purpose
!-------
!oxidation of SO2 and uptake of other gases in the aqueous phase
!
!interface
!---------
!PY0    initial concentration - at the end of this routine PY0 will be updated 
!       with the final concentrations 'PY' because otherwise the conc.changes due 
!       to wet chemistry will be 'forgotten' in TM5_DO_EBI. 
!PDT    chemistry timestep
!PY     concentrations at time is t
!   
!method
!------
!implicit solution of oxidation of SO2
!
!external
!--------
!none
!
!reference
!---------
!-
!**********************************************************************

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RNAVO,RMD
USE YOM_YGFL , ONLY : TYPE_GFLD
USE TM5_CHEM_MODULE,ONLY: NREACW,NTLOW,NTEMP,KSO2HP,KSO2O3, &
 &  HENRY,IMSA,ISO4,IHNO3,INH4,IH2O2,IO3,ISO2,INH3,INO3_A

IMPLICIT NONE

! input/output

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON
REAL(KIND=JPRB),INTENT(IN)         :: PDT 
REAL(KIND=JPRB), INTENT(INOUT)     :: PY0(KLON,YGFL%NCHEM)
REAL(KIND=JPRB), INTENT(IN)        :: PTEMP(KLON)
REAL(KIND=JPRB), INTENT(IN)        :: PRS(KLON)
REAL(KIND=JPRB), INTENT(IN)        :: PLP(KLON)
REAL(KIND=JPRB), INTENT(IN)        :: PAP(KLON) ! cloud cover
REAL(KIND=JPRB), INTENT(OUT)       :: PHPLUS(KLON) ! concentration H+
REAL(KIND=JPRB), INTENT(INOUT)     :: PY(KLON,YGFL%NCHEM+3) 


! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

INTEGER(KIND=JPIM)        :: ITEMP,ITER,JL
REAL(KIND=JPRB)           :: ZX1,ZX2,ZX3,ZDSO2,ZDNH3,ZXSO2O3A,ZXSO2O3B
REAL(KIND=JPRD)           :: ZX1_DP,ZX2_DP,ZX3_DP,ZSO2X_DP,ZDISC_DP,ZB1_DP,ZB2_DP
REAL(KIND=JPRB),PARAMETER :: ZCO2=3.90E-4 ! approx Global mean CO2 mixing ratio, 2000-2010 period
REAL(KIND=JPRB),PARAMETER :: ZRG=0.08314
REAL(KIND=JPRB)                 :: ZAIRD
REAL(KIND=JPRB),DIMENSION(KLON) :: ZXLIQ   ! liquid water content 
REAL(KIND=JPRB),DIMENSION(KLON) :: ZHKSO2  ! Henry's constant for sulfur dioxide
REAL(KIND=JPRB),DIMENSION(KLON) :: ZHKH2O2 ! Henry's constant for hydroperoxide
REAL(KIND=JPRB),DIMENSION(KLON) :: ZHKO3       ! Henry's constant for ozone
REAL(KIND=JPRB),DIMENSION(KLON) :: ZDKSO2      ! Dissociation constant for SO2
REAL(KIND=JPRB),DIMENSION(KLON) :: ZDKHSO3     ! Dissociation constant for HSO3-
REAL(KIND=JPRB),DIMENSION(KLON) :: ZDKH2O      ! dissociation constant water
REAL(KIND=JPRB),DIMENSION(KLON) :: ZDKNH3      ! dissociation constant ammonia
REAL(KIND=JPRB),DIMENSION(KLON) :: ZHKNH3      ! Henry's constant ammonia
REAL(KIND=JPRB),DIMENSION(KLON) :: ZHKCO2      ! Henry's constant CO2
REAL(KIND=JPRB),DIMENSION(KLON) :: ZDKCO2      ! Dissociation constant CO2
REAL(KIND=JPRB) :: ZPHS4                       ! effective dissolvation of S(IV)
REAL(KIND=JPRB) :: ZPHSO2                      ! effective dissolvation of SO2
REAL(KIND=JPRB) :: ZPHH2O2                     ! effective dissolvation of H2O2
REAL(KIND=JPRB) :: ZPHOZONE                    ! effective dissolvation of O3
REAL(KIND=JPRB),DIMENSION(KLON) :: ZPH         ! pH
REAL(KIND=JPRB) :: ZA1,ZA2,ZA,ZB,ZC,ZZ         ! help variables
REAL(KIND=JPRB) :: ZXCOV,ZXLIQR,ZXL,ZTEMP,ZRT,ZTR ! meteo
REAL(KIND=JPRB),DIMENSION(KLON,NREACW)   :: ZRW ! reaction rates
LOGICAL,DIMENSION(KLON)::  LLCLOUDY
REAL(KIND=JPRD),PARAMETER   :: Z_ONE=1._JPRD-1E-10_JPRD

! --- begin --------------------------------
IF (LHOOK) CALL DR_HOOK('TM5_WETCHEM',0,ZHOOK_HANDLE )

!-----------------------------
! wet phase reactions
!-----------------------------


DO JL=KIDIA,KFDIA
  LLCLOUDY(JL)=.FALSE.

   ZTEMP=PTEMP(JL)
   ZAIRD = 7.24291E16_JPRB*PRS(JL)/ZTEMP
   ! zx1: kg water to m3 water  (m3 water/ kg air)
   ZX1=PLP(JL)*1.E-3_JPRB 
   ZX2=RMD*1.E3_JPRB*ZAIRD/RNAVO    ! kg/m3 (air)
   ZXLIQ(JL) = ZX1/ZX2                            ! dimensionless number (m^3/m^3)
   ! avoid negatives and artificial values(1e-10 is about 0.0001 g/m3)
   IF ( ZXLIQ(JL) < 1E-10_JPRB ) ZXLIQ(JL)=0._JPRB 

  ! lwc is dimensionless 
  IF ((ZXLIQ(JL) > 1E-10_JPRB).AND.(PAP(JL) > 0.01_JPRB)) THEN 
    LLCLOUDY(JL)=.TRUE. 
    ZTR=(1._JPRB/ZTEMP-1._JPRB/298_JPRB)
    ZRT=ZTEMP*ZRG
    ITEMP=MIN(NTEMP,MAX(1_JPIM, NINT(ZTEMP-FLOAT(NTLOW))))
    !
    ! Henry and dissociation equilibria
    !
    ZDKH2O(JL) =1.01E-14_JPRB*EXP(-6706.0_JPRB *ZTR)   !h2o<=>hplus+so3--
    ZHKCO2(JL) =3.4E-2_JPRB*EXP(2420._JPRB*ZTR)        ! is already dimensionless
    !CMK BUG 092005 hkco2(JL) =3.4e-2*(2420.*ZTR)     ! is already dimensionless
    ZDKCO2(JL) =4.5E-7_JPRB*EXP(-1000._JPRB*ZTR)       !co2aq<=>hco3- + hplus
    ZHKSO2(JL) =HENRY(ISO2,ITEMP)*ZRT                   !dimensionless
    ZDKNH3(JL) =1.7E-5_JPRB*EXP(-450._JPRB*ZTR)        !nh3<=>nh4+ + OH-
    ZHKNH3(JL) =HENRY(INH3,ITEMP)*ZRT                   !dimensionless
    ZHKH2O2(JL)=HENRY(IH2O2,ITEMP)*ZRT                  !dimensionless
    ZHKO3(JL)  =HENRY(IO3,ITEMP)*ZRT                    !dimensionless
    ZDKSO2(JL) =1.3E-2_JPRB*EXP(2090._JPRB*ZTR)        !so2<=>hso3m+hplus
    ZDKHSO3(JL)=6.6E-8_JPRB*EXP(1510._JPRB*ZTR)        !hso3m<=>so3-- + hplus
    !
    ! calculate H+ from initial sulfate, ammonium, hno3, and nh3
    ! if solution is strongly acidic no further calculations are performed
    !

    ZXL=ZXLIQ(JL)*RNAVO*1E-3_JPRB/PAP(JL)
    !ZX1 is initial strong acidity from SO4 and NO3
    !
    !acidity from strong acids alone
    !  
    PHPLUS(JL)=(2._JPRB*PY0(JL,ISO4)+PY0(JL,IMSA)-PY0(JL,INH4)+&
       &               PY0(JL,IHNO3)+PY0(JL,INO3_A))/ZXL
    PHPLUS(JL)=MAX(PHPLUS(JL),1E-20_JPRB)
 
  ENDIF
ENDDO
DO ITER=1,10
  DO JL=KIDIA,KFDIA
         ! only if solution pH>4.5
       IF ( LLCLOUDY(JL)  ) THEN
          IF (  (PHPLUS(JL) < 3E-5_JPRB) ) THEN
            ZXL=ZXLIQ(JL)*RNAVO*1E-3_JPRB/PAP(JL)
            ZX1=(2._JPRB*PY0(JL,ISO4)+PY0(JL,IMSA)+PY0(JL,IHNO3)+PY0(JL,INO3_A))/ZXL    
            !ZX2 is initial total NHx
            ZX2=(PY0(JL,INH3)+PY0(JL,INH4))/ZXL
            !ZX3 is combined dissolution and solubility const for CO2
            ZX3=ZDKCO2(JL)*ZHKCO2(JL)*ZCO2 
            ZA1=ZDKH2O(JL)/ZDKNH3(JL)*(1._JPRB+1._JPRB/ZHKNH3(JL)) ! integration constant
            ZA2=PY0(JL,ISO2)/ZXL      !initial SO2
            ZZ =ZA2/(PHPLUS(JL)/ZDKSO2(JL)*(1._JPRB+1._JPRB/ZHKSO2(JL))+&
              &  ZDKHSO3(JL)/PHPLUS(JL)+1._JPRB)
            ZA=1._JPRB+ZX2/(ZA1+PHPLUS(JL))
            ZB=-ZX1-ZZ
            ZC=-ZX3-2._JPRB*ZDKHSO3(JL)*ZZ
            ZZ=MAX(0._JPRB,(ZB*ZB-4._JPRB*ZA*ZC))
            PHPLUS(JL)=MAX(1.E-10_JPRB,(-ZB+SQRT(ZZ))/(2._JPRB*ZA))
        ENDIF
      ENDIF
   ENDDO ! 
ENDDO  !iter
DO JL=KIDIA,KFDIA
  IF (LLCLOUDY(JL)) THEN
    ZTEMP=PTEMP(JL)
    ZTR=(1._JPRB/ZTEMP-1._JPRB/298_JPRB)
    ZXLIQR=ZXLIQ(JL)/PAP(JL)
    ZXL=ZXLIQ(JL)*RNAVO*1E-3_JPRB/PAP(JL)
    ZPH(JL)=-LOG10(PHPLUS(JL))     ! pH for diagnostics 

    ! phase factor ratio of aqueous phase to gas phase concentration

    ZPHS4   = ZHKSO2(JL) *(1._JPRB+ZDKSO2(JL)/PHPLUS(JL)+&
            & ZDKSO2(JL)*ZDKHSO3(JL)/PHPLUS(JL)/PHPLUS(JL))*ZXLIQR
    ZPHSO2  =ZHKSO2(JL)  *ZXLIQR
    ZPHH2O2 =ZHKH2O2(JL) *ZXLIQR
    ZPHOZONE=ZHKO3(JL)   *ZXLIQR

    ! the original rate equations could be partly in look-up table

    ZRW(JL,KSO2HP) =8E4_JPRB    *EXP(-3560._JPRB*ZTR)/(0.1_JPRB+PHPLUS(JL))
    ZXSO2O3A       =4.39E11_JPRB*EXP(-4131._JPRB/ZTEMP)+2.56E3_JPRB*EXP(-966._JPRB/ZTEMP)  !S(IV)
    ZXSO2O3B       =2.56E3_JPRB *EXP(-966._JPRB /ZTEMP)/PHPLUS(JL)  
    !divide by [H+]!S(IV)

    !  make rate constants dimensionless by multiplying 
    !  by (1./ZXLIQR/avo=6e20)  
    !  multiply with the fractions of concentrations residing 
    !  in the aqueous phase

    ZRW(JL,KSO2HP)=ZRW(JL,KSO2HP)/ZXL*ZPHSO2/(1._JPRB+ZPHS4)*ZPHH2O2/(1._JPRB+ZPHH2O2)
    ZRW(JL,KSO2O3)=(ZXSO2O3A+ZXSO2O3B)/ZXL*ZPHS4/(1._JPRB+ZPHS4)*ZPHOZONE/&
                & (1._JPRB+ZPHOZONE)
  ENDIF !cloudy
ENDDO ! JL LOOP


DO JL=KIDIA,KFDIA
   !
   ! only cloud chemistry if substantial amount of clouds are present 
   !
   IF (LLCLOUDY(JL)) THEN
      !
      ! oxidation of S(IV) by O3
      !
      ZSO2X_DP=PY0(JL,ISO2)
      ZXCOV=PAP(JL)
      ZX1=MIN(100._JPRB,ZRW(JL,KSO2O3)*PY0(JL,IO3)*PDT)
      ZDSO2=PY0(JL,ISO2)*ZXCOV*(EXP(-ZX1)-1._JPRB)
      !only applied to ZXCOV part of cloud
      ZDSO2=MAX(-PY0(JL,IO3)*ZXCOV,ZDSO2)! limit to o3 availability
      PY(JL,ISO2)=PY0(JL,ISO2)+ZDSO2 
      !NOTE CMK: paralel MPI should take care here!
      PY(JL,ISO4)=PY0(JL,ISO4)-ZDSO2
      PY(JL,IO3) =PY0(JL,IO3) +ZDSO2

      !
      ! oxidation of S(IV) by H2O2
      !
      !*** here we explicitly solve the dv: 
      !    y'=P-Q*y-R*y*y (P and Q are 0=>b3=0.)
      !
      ZSO2X_DP=PY(JL,ISO2)
      IF ( ZSO2X_DP > 1E-22_JPRD ) THEN
        ZB1_DP=MAX(1E-22_JPRD,ZRW(JL,KSO2HP))
        ZB2_DP=ZB1_DP*(PY0(JL,IH2O2)-ZSO2X_DP)
        ZDISC_DP=MIN(1._JPRB,SQRT(ZB2_DP*ZB2_DP))        ! disc is b2 for b3=0.0
        ZX1_DP=(ZB2_DP-ZDISC_DP)/(-2._JPRB*ZB1_DP)          ! in this case ZX1 =0.
        ZX2_DP=(ZB2_DP+ZDISC_DP)/(-2._JPRB*ZB1_DP)
        ZX3_DP=MIN(Z_ONE,(ZSO2X_DP-ZX1_DP)/(ZSO2X_DP-ZX2_DP)*EXP(-ZDISC_DP*PDT))
        ZSO2X_DP=(ZX1_DP-ZX2_DP*ZX3_DP)/(1._JPRD-ZX3_DP)
        ZDSO2=(ZSO2X_DP-PY(JL,ISO2))*ZXCOV
        ZDSO2=MAX(ZDSO2,-PY0(JL,IH2O2)*ZXCOV)
        PY(JL,ISO2)  = PY(JL,ISO2)+ZDSO2         ! dso2 is loss of SO2 and H2O2
        PY(JL,ISO4)  = PY(JL,ISO4)-ZDSO2
        PY(JL,IH2O2) = PY0(JL,IH2O2)+ZDSO2 
      ENDIF

      !
      ! NH3 uptake in cloud droplets is limited by H2SO4 availability
      ! no HNO3 is considered at this point
      ! assume instantaneous uptake of NH3 incloud  only in cloudy part
      !
      ZDNH3=MAX((2._JPRB*PY(JL,ISO4)+PY0(JL,IMSA)-PY0(JL,INH4))*ZXCOV,0._JPRB)
      ZDNH3=MAX(-PY0(JL,INH3)*ZXCOV,-ZDNH3)
      PY(JL,INH3)=PY0(JL,INH3)+ZDNH3    ! dnh3 is loss of NH3  
      PY(JL,INH4)=PY0(JL,INH4)-ZDNH3

      ! Finally update PY0 concentrations
     
      
      ! Finally also update PY0 concentrations
      
      PY0(JL,INH3)=PY(JL,INH3)
      PY0(JL,INH4)=PY(JL,INH4)
      PY0(JL,ISO2)=PY(JL,ISO2)
      PY0(JL,ISO4)=PY(JL,ISO4)
      PY0(JL,IO3) =PY(JL,IO3)
      PY0(JL,IH2O2)=PY(JL,IH2O2)

   ENDIF     !cloudy
ENDDO 


IF (LHOOK) CALL DR_HOOK('TM5_WETCHEM',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_WETCHEM
