! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_WETCHEM_POINT(YGFL,PDT,PTEMP,PAP,PRS,PLP,PHPLUS,PY)
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

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RNAVO,RMD
USE YOM_YGFL , ONLY : TYPE_GFLD
USE TM5_CHEM_MODULE,ONLY: NREACW,NTLOW,NTEMP,KSO2HP,KSO2O3
USE BASCOETM5_MODULE, ONLY : HENRY, & 
        & IMSA,ISO4,IHNO3,INH4,IH2O2,IO3,ISO2,INH3,INO3_A

IMPLICIT NONE

! input/output

TYPE(TYPE_GFLD),INTENT(INOUT)      :: YGFL
REAL(KIND=JPRB),INTENT(IN)         :: PDT 
REAL(KIND=JPRB), INTENT(IN)        :: PTEMP
REAL(KIND=JPRB), INTENT(IN)        :: PRS
REAL(KIND=JPRB), INTENT(IN)        :: PLP
REAL(KIND=JPRB), INTENT(IN)        :: PAP ! cloud cover
REAL(KIND=JPRB), INTENT(OUT)       :: PHPLUS ! concentration H+
REAL(KIND=JPRB), INTENT(INOUT)     :: PY(YGFL%NCHEM) 


! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

INTEGER(KIND=JPIM)        :: ITEMP,ITER
REAL(KIND=JPRB)           :: ZX1,ZX2,ZX3,ZB1,ZB2,ZSO2X,ZDSO2,ZDISC,ZDNH3,ZXSO2O3A,ZXSO2O3B
REAL(KIND=JPRB),PARAMETER :: ZCO2=3.90E-4 ! approx global mean CO2 mixing ratio, 2000-2010 period
REAL(KIND=JPRB),PARAMETER :: ZRG=0.08314
REAL(KIND=JPRB)           :: ZAIRD
REAL(KIND=JPRB)           :: ZXLIQ   ! liquid water content 
REAL(KIND=JPRB)           :: ZHKSO2  ! Henry's constant for sulfur dioxide
REAL(KIND=JPRB)           :: ZHKH2O2 ! Henry's constant for hydroperoxide
REAL(KIND=JPRB)           :: ZHKO3   ! Henry's constant for ozone
REAL(KIND=JPRB)           :: ZDKSO2  ! Dissociation constant for SO2
REAL(KIND=JPRB)           :: ZDKHSO3 ! Dissociation constant for HSO3-
REAL(KIND=JPRB)           :: ZDKH2O  ! dissociation constant water
REAL(KIND=JPRB)           :: ZDKNH3  ! dissociation constant ammonia
REAL(KIND=JPRB)           :: ZHKNH3  ! Henry's constant ammonia
REAL(KIND=JPRB)           :: ZHKCO2  ! Henry's constant CO2
REAL(KIND=JPRB)           :: ZDKCO2  ! Dissociation constant CO2
REAL(KIND=JPRB)           :: ZPHS4     ! effective dissolvation of S(IV)
REAL(KIND=JPRB)           :: ZPHSO2    ! effective dissolvation of SO2
REAL(KIND=JPRB)           :: ZPHH2O2   ! effective dissolvation of H2O2
REAL(KIND=JPRB)           :: ZPHOZONE  ! effective dissolvation of O3
REAL(KIND=JPRB)           :: ZPH       ! pH
REAL(KIND=JPRB)           :: ZA1,ZA2,ZA,ZB,ZC,ZZ         ! help variables
REAL(KIND=JPRB)           :: ZXCOV,ZXLIQR,ZXL,ZTEMP,ZRT,ZTR ! meteo
REAL(KIND=JPRB),DIMENSION(NREACW)   :: ZRW ! reaction rates
LOGICAL                   ::  LLCLOUDY

! --- begin --------------------------------
IF (LHOOK) CALL DR_HOOK('TM5_WETCHEM_POINT',0,ZHOOK_HANDLE )

!-----------------------------
! wet phase reactions
!-----------------------------


  LLCLOUDY=.FALSE.

  ZTEMP=PTEMP
  ZAIRD = 7.24291E16_JPRB*PRS/ZTEMP
  ! zx1: kg water to m3 water  (m3 water/ kg air)
  ZX1=PLP*1.E-3_JPRB 
  ZX2=RMD*1.E3_JPRB*ZAIRD/RNAVO    ! kg/m3 (air)
  ZXLIQ = ZX1/ZX2                            ! dimensionless number (m^3/m^3)
  ! avoid negatives and artificial values(1e-10 is about 0.0001 g/m3)
  IF ( ZXLIQ < 1E-10_JPRB ) ZXLIQ=0._JPRB 

  ! lwc is dimensionless 
  IF ((ZXLIQ > 1E-10_JPRB).AND.(PAP > 0.01_JPRB)) THEN 
    ! LLCLOUDY=.TRUE. 
    !
    ! calculate H+ from initial sulfate, ammonium, hno3, and nh3
    ! if solution is strongly acidic no further calculations are performed
    !

    ZXL=ZXLIQ*RNAVO*1E-3_JPRB/PAP
    !ZX1 is initial strong acidity from SO4 and NO3
    !
    !acidity from strong acids alone
    !  
    PHPLUS=(2._JPRB*PY(ISO4)+PY(IMSA)-PY(INH4)+&
       &               PY(IHNO3)+PY(INO3_A))/ZXL
    PHPLUS=MAX(PHPLUS,1E-20_JPRB)

    !VH directly Continue:  ENDIF
    !VH directly continue:  IF ( LLCLOUDY  ) THEN

    ZTR=(1._JPRB/ZTEMP-1._JPRB/298_JPRB)
    ZRT=ZTEMP*ZRG
    ITEMP=MIN(NTEMP,MAX(1_JPIM, NINT(ZTEMP-FLOAT(NTLOW))))
    !
    ! Henry and dissociation equilibria
    !
    ZDKH2O =1.01E-14_JPRB*EXP(-6706.0_JPRB *ZTR)   !h2o<=>hplus+so3--
    ZHKCO2 =3.4E-2_JPRB*EXP(2420._JPRB*ZTR)        ! is already dimensionless
    !CMK BUG 092005 hkco2 =3.4e-2*(2420.*ZTR)     ! is already dimensionless
    ZDKCO2 =4.5E-7_JPRB*EXP(-1000._JPRB*ZTR)       !co2aq<=>hco3- + hplus
    ZHKSO2 =HENRY(ISO2,ITEMP)*ZRT                   !dimensionless
    ZDKNH3 =1.7E-5_JPRB*EXP(-450._JPRB*ZTR)        !nh3<=>nh4+ + OH-
    ZHKNH3 =HENRY(INH3,ITEMP)*ZRT                   !dimensionless
    ZHKH2O2=HENRY(IH2O2,ITEMP)*ZRT                  !dimensionless
    ZHKO3  =HENRY(IO3,ITEMP)*ZRT                    !dimensionless
    ZDKSO2 =1.3E-2_JPRB*EXP(2090._JPRB*ZTR)        !so2<=>hso3m+hplus
    ZDKHSO3=6.6E-8_JPRB*EXP(1510._JPRB*ZTR)        !hso3m<=>so3-- + hplus

    ZXL=ZXLIQ*RNAVO*1E-3_JPRB/PAP
    ZX1=(2._JPRB*PY(ISO4)+PY(IMSA)+PY(IHNO3)+PY(INO3_A))/ZXL    
    !ZX2 is initial total NHx
    ZX2=(PY(INH3)+PY(INH4))/ZXL
    !ZX3 is combined dissolution and solubility const for CO2
    ZX3=ZDKCO2*ZHKCO2*ZCO2 
    ZA1=ZDKH2O/ZDKNH3*(1._JPRB+1._JPRB/ZHKNH3) ! integration constant
    ZA2=PY(ISO2)/ZXL      !initial SO2

    DO ITER=1,10
      IF (  (PHPLUS < 3E-5_JPRB) ) THEN
        ZZ =ZA2/(PHPLUS/ZDKSO2*(1._JPRB+1._JPRB/ZHKSO2)+&
          &  ZDKHSO3/PHPLUS+1._JPRB)
        ZA=1._JPRB+ZX2/(ZA1+PHPLUS)
        ZB=-ZX1-ZZ
        ZC=-ZX3-2._JPRB*ZDKHSO3*ZZ
        ZZ=MAX(0._JPRB,(ZB*ZB-4._JPRB*ZA*ZC))
        PHPLUS=MAX(1.E-10_JPRB,(-ZB+SQRT(ZZ))/(2._JPRB*ZA))
      ENDIF
    ENDDO  !iter

    ZXLIQR=ZXLIQ/PAP
    ZXL=ZXLIQ*RNAVO*1E-3_JPRB/PAP
    ZPH=-LOG10(PHPLUS)     ! pH for diagnostics 

    ! phase factor ratio of aqueous phase to gas phase concentration

    ZPHS4   = ZHKSO2 *(1._JPRB+ZDKSO2/PHPLUS+&
            & ZDKSO2*ZDKHSO3/PHPLUS/PHPLUS)*ZXLIQR
    ZPHSO2  =ZHKSO2  *ZXLIQR
    ZPHH2O2 =ZHKH2O2 *ZXLIQR
    ZPHOZONE=ZHKO3   *ZXLIQR

    ! the original rate equations could be partly in look-up table

    ZRW(KSO2HP) =8E4_JPRB    *EXP(-3560._JPRB*ZTR)/(0.1_JPRB+PHPLUS)
    ZXSO2O3A       =4.39E11_JPRB*EXP(-4131._JPRB/ZTEMP)+2.56E3_JPRB*EXP(-966._JPRB/ZTEMP)  !S(IV)
    ZXSO2O3B       =2.56E3_JPRB *EXP(-966._JPRB /ZTEMP)/PHPLUS  
    !divide by [H+]!S(IV)

    !  make rate constants dimensionless by multiplying 
    !  by (1./ZXLIQR/avo=6e20)  
    !  multiply with the fractions of concentrations residing 
    !  in the aqueous phase

    ZRW(KSO2HP)=ZRW(KSO2HP)/ZXL*ZPHSO2/(1._JPRB+ZPHS4)*ZPHH2O2/(1._JPRB+ZPHH2O2)
    ZRW(KSO2O3)=(ZXSO2O3A+ZXSO2O3B)/ZXL*ZPHS4/(1._JPRB+ZPHS4)*ZPHOZONE/&
                & (1._JPRB+ZPHOZONE)

   !
   ! only cloud chemistry if substantial amount of clouds are present 
   !
   !
   ! oxidation of S(IV) by O3
   !
   ZSO2X=PY(ISO2)
   ZXCOV=PAP
   ZX1=MIN(100._JPRB,ZRW(KSO2O3)*PY(IO3)*PDT)
   ZDSO2=PY(ISO2)*ZXCOV*(EXP(-ZX1)-1._JPRB)
   !only applied to ZXCOV part of cloud
   ZDSO2=MAX(-PY(IO3)*ZXCOV,ZDSO2)! limit to o3 availability
   PY(ISO2)=PY(ISO2)+ZDSO2 
   !NOTE CMK: paralel MPI should take care here!
   PY(ISO4)=PY(ISO4)-ZDSO2
   PY(IO3) =PY(IO3) +ZDSO2

   !
   ! oxidation of S(IV) by H2O2
   !
   !*** here we explicitly solve the dv: 
   !    y'=P-Q*y-R*y*y (P and Q are 0=>b3=0.)
   !
   ZSO2X=PY(ISO2)
   IF ( ZSO2X > TINY(ZSO2X) ) THEN
     ZB1=ZRW(KSO2HP)
     ZB2=ZB1*(PY(IH2O2)-ZSO2X)
     ZDISC=MIN(100._JPRB,SQRT(ZB2*ZB2))     ! disc is b2 for b3=0.0
     ZX1=(ZB2-ZDISC)/(-2._JPRB*ZB1)         ! in this case ZX1 =0.
     ZX2=(ZB2+ZDISC)/(-2._JPRB*ZB1)
     ZX3=(ZSO2X-ZX1)/(ZSO2X-ZX2)*EXP(-ZDISC*PDT)
     IF (ZX3 == 1.0_JPRB ) ZX3=0.999_JPRB
     ZSO2X=(ZX1-ZX2*ZX3)/(1.-ZX3)
     ZDSO2=(ZSO2X-PY(ISO2))*ZXCOV
     ZDSO2=MAX(ZDSO2,-PY(IH2O2)*ZXCOV)
     PY(ISO2)  = PY(ISO2)+ZDSO2             ! dso2 is loss of SO2 and H2O2
     PY(ISO4)  = PY(ISO4)-ZDSO2
     PY(IH2O2) = PY(IH2O2)+ZDSO2 
   ENDIF

   !
   ! NH3 uptake in cloud droplets is limited by H2SO4 availability
   ! no HNO3 is considered at this point
   ! assume instantaneous uptake of NH3 incloud  only in cloudy part
   !
   ZDNH3=MAX((2._JPRB*PY(ISO4)+PY(IMSA)-PY(INH4))*ZXCOV,0._JPRB)
   ZDNH3=MAX(-PY(INH3)*ZXCOV,-ZDNH3)
   PY(INH3)=PY(INH3)+ZDNH3    ! dnh3 is loss of NH3  
   PY(INH4)=PY(INH4)-ZDNH3

   
   ! Finally also update PY concentrations
   
   ! PY0(INH3)=PY(INH3)
   ! PY0(INH4)=PY(INH4)
   ! PY0(ISO2)=PY(ISO2)
   ! PY0(ISO4)=PY(ISO4)
   ! PY0(IO3) =PY(IO3)
   ! PY0(IH2O2)=PY(IH2O2)

ENDIF     !cloudy


IF (LHOOK) CALL DR_HOOK('TM5_WETCHEM_POINT',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_WETCHEM_POINT
