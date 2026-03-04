! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AERC_SCAVL19 &
 & ( YDMODEL,KIDIA , KFDIA  , KLON , KLEV , KAER , KAER_SCAV, KSTEP, PTSPHY,&
 &   PRSF1, PDP , PTP, PFLXR, PFLXS, PCLCOV, PCLWAT, PCLICE, PRAIN, PSNOW, PCEN, PTENC0, &
 &   PTENC1, PFAERO, PPRCOV  )

!*** * AERC_SCAVL19* - IN-CLOUD AND BELOW CLOUD SCAVENGING OF TRACERS
!      CALLED SEPARATELY FOR CONVECTIVE AND LARGE-SCALE PRECIP. 
! COPY of CHEM_SCAV 

! INPUTS:
! -------
! KSTEP :  Time step number
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KLEV  :  Number of Levels
! KCHEM :  Number of chemistry tracers 
! KCHEM :  Number of chemistry tracers with wet depostion


! PTSPHY:  Time step length in seconds
! PDP(KLON,KLEV)              :  PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRSF1(KLON,KLEV)            :  Mid-level pressure           (Pa)
! PTP(KLON,KLEV)              :  Temperature in            (T)
! PCLWAT  (KLON,KLEV)         :  Cloud water content    (kg/kg) - stratiform and convective 
! PCLICE  (KLON,KLEV)         :  Cloud ice water content    (kg/kg) - stratiform and convective
! PRAIN (KLON,KLEV)           :  Rain water content    (kg/kg) for stratiform precip
! PSNOW (KLON,KLEV)           :  Snow water content    (kg/kg) for stratiform precip
! PFLXR   (KLON,KLEV+1)       :  Precip flux rain       (kg/m2s) - either stratiform or convective (depending on argument)
! PFLXS   (KLON,KLEV+1)       :  Precip flux snow       (kg/m2s) - either stratoform or convective (depending on argument)
! PCLCOV  (KLON,KLEV)         :  cloud fraction   0..1
! PPRCOV  (KLON,KLEV)         :  precipitation fraction   0..1
! PCEN(KLON,KLEV,KCHEM)       :  CONCENTRATION OF TRACERS           (kg/kg)
! PTENC0(KLON,KLEV,KCHEM)     :  TOTAL TENDENCY OF CONCENTRATION OF TRACERS BEFORE(kg/kg s-1)
!
! NB: PCLWAT is the in-cloud water mixing ratio
! OUTPUTS:
! -------
! PTENC1 (KLON,KLEV,KCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS after (kg/kg s-1)
!
!**   INTERFACE.
!     ----------
!          *AERC_SCAVL19* IS CALLED FROM *AER_PHY3*.
!
!     AUTHOR.
!     -------
!        Johannes Flemming 
!        
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2009-11-09
!        J. Flemming : use precip and cloud cover area: 2012-10-10 
!        J. Flemming : check if trandfer or henry limited rain out (see Jacobs 2000)  2013-6-10 
!        J. Flemming : use stratifrom rain and snow water content: 2013-6-13 
!-----------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RG, RD, RLVTT ,RLSTT ,RTT 
USE YOETHF   , ONLY :  R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,&
 & RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R


IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

TYPE(MODEL)       ,INTENT(INOUT):: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KAER, KAER_SCAV, KSTEP

REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV) , PRSF1(KLON,KLEV) , PTP(KLON,KLEV)   
REAL(KIND=JPRB),INTENT(IN)    :: PCLCOV(KLON,KLEV) , PCLWAT(KLON,KLEV),PCLICE(KLON,KLEV) , PRAIN(KLON,KLEV), PSNOW(KLON,KLEV)  
REAL(KIND=JPRB),OPTIONAL, INTENT(IN)    ::  PPRCOV(KLON,KLEV)  
REAL(KIND=JPRB),INTENT(IN)    :: PFLXR(KLON,KLEV+1), PFLXS(KLON,KLEV+1)
REAL(KIND=JPRB),INTENT(IN)    :: PTENC0(KLON,KLEV,KAER), PCEN(KLON,KLEV,KAER)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY

REAL(KIND=JPRB),INTENT(OUT)   :: PTENC1(KLON,KLEV,KAER)
REAL(KIND=JPRB),INTENT(OUT)   :: PFAERO(KLON,KAER)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM), PARAMETER :: IMODE = 1 ! Choice of parameter scheme. Scheme 2 is optimized for TM5
INTEGER(KIND=JPIM) :: JK, JL, JT, IWETDEP

REAL(KIND=JPRB) :: ZSCAV
REAL(KIND=JPRB) :: ZBETA,  ZBETAR, ZBETARI, ZBETASI
REAL(KIND=JPRB) :: ZDP, ZDX, ZDZ, ZFRAC, ZRHO, ZFUNC
REAL(KIND=JPRB) :: ZEPSQLIQ, ZEPSFLX, ZEPSQLIQ2
REAL(KIND=JPRB) :: ZRD, ZLMMR2VMR, ZRET, ZHNRYEFT
REAL(KIND=JPRB) :: ZLIQ2TOT, ZICE2TOT, ZLIQ2GAS, ZICE2GAS, ZRAINW, ZFALLSP
REAL(KIND=JPRB) :: ZPRCOV, ZCLCOV, ZFLXR, ZFLXS,ZFLXRB, ZFLXSB,  ZCLWAT, ZCLICE,ZCLTOT, ZMINCLCOV, ZMINPRCOV 
REAL(KIND=JPRB) :: ZCLTOTW, ZCLTOTI
 ! Interstitial Fraction: 30% of aerosol remains in atmosphere
REAL(KIND=JPRB)  ::  ZINTERST_FR  
REAL(KIND=JPRB),PARAMETER  :: ZDGHNO3 = 0.136      ! viscosity of HNO3 in [cm2/s] 
REAL(KIND=JPRB),PARAMETER  :: ZDGAIR  = 0.133      ! viscosity of air in [cm2/s] 
REAL(KIND=JPRB),PARAMETER  :: ZXMHNO3    =1.008_JPRB + 14.007_JPRB + 3*16.0_JPRB ! HNO3 tracer mass
REAL(KIND=JPRB),PARAMETER :: ZHPLUS =3.16227E-6_JPRB ! is rain water ph=5.5 H+ 
REAL(KIND=JPRB)  :: ZRLWC,ZRDRAD, ZRU, ZNRE,ZNSC,ZNSH,  ZKG, ZTR, ZKSO2, ZKHSO3, ZFACTSO2
REAL(KIND=JPRB)            :: ZRL      ! composite factor of Rgas and liquid water content of raining cloud
                                       ! rgas (8.314 J/mol/K) ---> 0.08314 atm/(mol/l)/K
                                       ! 1e-6 corresponds to 1 g/m3 dimensionless 
REAL(KIND=JPRB)   :: ZWLOSS(KLON,KAER)
REAL(KIND=JPRB)   :: ZKMIN, ZKNEW, ZF, ZKMINICE

LOGICAL :: LLWDAER 
LOGICAL :: LLPRINT, LLCHEM_WDFR, LLCONV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttre.func.h"
!#include "fccld.h"


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AERC_SCAVL19',0,ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP,YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM, &
 & YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM)
ASSOCIATE(YCHEM=>YGFL%YCHEM,YAERO_DESC=>YDEAERATM%YAERO_DESC, &
 & LCHEM_WDFR=>YDCHEM%LCHEM_WDFR, &
 & RCOVPMIN=>YDECLDP%RCOVPMIN, &
 & RFRAER=>YDEAERSNK%RFRAER)
LLPRINT=.FALSE.
LLCHEM_WDFR=LCHEM_WDFR
LLCONV=.FALSE.

ZINTERST_FR = 0.3_JPRB
ZKMIN=1.E-4_JPRB
ZKMINICE=1.E-6_JPRB
IF ( LLCHEM_WDFR ) ZINTERST_FR = 0.0_JPRB 
  
! deposition for convective Precip
IF (.NOT. PRESENT(PPRCOV) )  LLCONV=.TRUE.

!* set flux to zero 
ZWLOSS(:,:)=0._JPRB

ZEPSFLX =1.E-18_JPRB
ZEPSQLIQ=2.E-6_JPRB
ZEPSQLIQ2=2.E-7_JPRB

!* mini cloud cover and precip cover
ZMINCLCOV=0.001_JPRB
ZMINPRCOV=0.001_JPRB

!* precip fall speed 
ZFALLSP = 5.0_JPRB  

!* ideal gas constant in atm M-1
ZRD=1000.0_JPRB * RD * 9.8692_JPRB / 1000000.0_JPRB
ZRD=0.082_JPRB
ZRL=RG/1E2_JPRB*1E-6_JPRB 
PFAERO(:,:)=0._JPRB

!* initialisation    
!* re-evaporation fraction to account for drop shrinking without releasing species
! Jacob says 
!ZFRAC=0.5_JPRB
ZFRAC=0.2_JPRB

!* update tendecies  

!* if not LLCHEM_WDFR 
ZPRCOV =  1.0_JPRB
ZCLCOV =  1.0_JPRB

PTENC1(KIDIA:KFDIA,1:KLEV,1:KAER) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KAER)

!* LOOP OVER LAYERS  
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

    !- precip flux at box top
    ZFLXR =  PFLXR(JL,JK)
    ZFLXS =  PFLXS(JL,JK)
    ZFLXRB =  PFLXR(JL,JK+1)
    ZFLXSB =  PFLXS(JL,JK+1)

!--  Precip flux change grid box average - indicator for rain/snow formation or
!evaporation
    ZBETA = ( ZFLXRB - ZFLXR ) + ( ZFLXSB - ZFLXS  )  ! in kg m-2 s-1
    ZBETARI = ZFLXRB - ZFLXR
    ZBETASI = ZFLXSB - ZFLXS
    ZCLWAT =  PCLWAT(JL,JK)
    ZCLICE =  PCLICE(JL,JK)

    IF (  ZBETA ==  0.0_JPRB .AND. ZFLXR < ZEPSFLX  )  CYCLE
!     IF ( LLCHEM_WDFR ) THEN
!* use effective flux and cloud cover for liquid water content
      IF (.NOT. LLCONV ) THEN
        ZPRCOV = MIN(MAX(RCOVPMIN,PPRCOV(JL,JK)),1.0_JPRB)  ! grid scale precip
        ZCLCOV = MIN(1.0_JPRB,MAX(0.0_JPRB,PCLCOV(JL,JK)))
      ELSE
        ZPRCOV = 0.05_JPRB    ! convective precip
        ZCLCOV = 0.05_JPRB    ! convective precip
      ENDIF
      IF (ZPRCOV >= ZMINPRCOV) THEN
        ZFLXR =  ZFLXR / ZPRCOV
        ZFLXS =  ZFLXS / ZPRCOV
        ZFLXRB =  ZFLXRB / ZPRCOV
       ZFLXSB =  ZFLXSB / ZPRCOV
        ZBETA = ZBETA / ZPRCOV
      ELSE
        ZFLXR =  0.0_JPRB
        ZFLXS =  0.0_JPRB
        ZBETA =  0.0_JPRB
      ENDIF
      IF (ZCLCOV >= ZMINCLCOV) THEN
        ZCLWAT =  MAX(0.0_JPRB,PCLWAT(JL,JK) )
        ZCLICE =  MAX(0.0_JPRB,PCLICE(JL,JK) )
      ELSE
        ZCLWAT = 0.0_JPRB
        ZCLICE = 0.0_JPRB
      ENDIF
!    ENDIF

!* calculate air density, layer depth 
     ZRHO=PRSF1(JL,JK)/(RD*PTP(JL,JK))
     ZDZ= PDP(JL,JK) / (ZRHO*RG)
     ZDP =  PDP(JL,JK)
     ZBETARI=ZBETARI*RG/ZDP  ! in s-1
     ZBETASI=ZBETASI*RG/ZDP  ! in s-1
! clwc mass mixing ratio to volume mixing ratio : rho air / rho cloud water - note it is liquid water not vapor
    ZLMMR2VMR=ZRHO/1000.0_JPRB 
    IWETDEP=0

!* LOOP over species
    DO JT=1,KAER

       LLWDAER=.TRUE. 
       IWETDEP=IWETDEP+1

       IF ( YAERO_DESC(JT)%CNAME == 'Sulphate_SO2'  ) THEN
          LLWDAER = .FALSE.
!* Retention coefficient for mixed clouds
         IF (PTP(JL,JK) < 268.0_JPRB) THEN
            ZRET=0.02_JPRB   ! For others assume retention of 0.02
         ELSE
!* otherwise 1.0  
           ZRET=1.0_JPRB
         ENDIF

!* Henry coefficicnet 
         ! Henrys are in M/atm = (mol/liter)/atm - RD has to be changed 
         ! Temperature dependency of henry constant 
         ZTR=(1.0_JPRB/PTP(JL,JK)-1.0_JPRB/298.0_JPRB)
!! effective Henry coeff for SO2 (Seinfeld Pandis, 1998. p 350)
         ZKSO2=1.7E-2_JPRB*EXP(2090.0_JPRB*ZTR)        ! so2<=>hso3m+hplus
         ZKHSO3=6.6E-8_JPRB*EXP(1510.0_JPRB*ZTR)     ! hso3m<=>so3-- + hplus
         ZFACTSO2=1.0_JPRB + ZKSO2 / ZHPLUS  + (ZKSO2*ZKHSO3)/(ZHPLUS**2.0_JPRB)
         ZHNRYEFT=ZFACTSO2 * 1.23_JPRB * EXP ( 3000.0_JPRB * ZTR)

!* liquid phase to gas phase fraction - henry's law
         IF (IMODE == 1_JPIM) THEN
         !* Standard - (this leads to too low scavenging rates in comparison to TM5)
           ZLIQ2GAS=ZHNRYEFT*ZCLWAT*ZRD*PTP(JL,JK)*ZLMMR2VMR  ! ZCLWAT in kg/kg  
         ELSEIF(IMODE == 2_JPIM) THEN
           !*VH - compute liq/gas fraction on old-fashioned TM5 way, 
           !*VH - using constant liquid water content, part of ZRL
           ZLIQ2GAS=ZHNRYEFT*PTP(JL,JK)*ZRL
         ENDIF
!* ice phase to gas phase fraction set to zero appart from HNO3 and H2O2
         ZICE2GAS=0.0_JPRB 
!* liquid phase to total fraction 
         ZLIQ2TOT=ZLIQ2GAS/(1.0_JPRB+ZLIQ2GAS+ZICE2GAS)  
!* ice phase to total fraction 
         ZICE2TOT=ZICE2GAS/(1.0_JPRB+ZLIQ2GAS+ZICE2GAS)  
       ENDIF ! SO2 

!* total cloud water ice
        ZCLTOTW=ZCLWAT + PTSPHY * ZBETARI
        ZCLTOTI=ZCLICE + PTSPHY * ZBETASI
        ZCLTOT=ZCLTOTW+ZCLTOTI
        ZSCAV=0._JPRB
        ZFUNC=0._JPRB
!* Rain-out in Cloud
        IF (ZCLTOTW > ZEPSQLIQ .AND. ZBETARI > 0._JPRB .AND. YAERO_DESC(JT)%RSCAVIN > 0.0_JPRB ) THEN
!* relative water loss due to precip formation

        !* relative water loss due to precipA  formation
          ZKNEW = ZKMIN + ZBETARI/ZCLTOTW
          ZF=ZCLCOV * ZBETARI/ (ZKNEW * ZCLTOTW)
          ZSCAV=(YAERO_DESC(JT)%RSCAVIN)*ZKNEW
          ZFUNC=(EXP(-ZSCAV*PTSPHY)-1._JPRB)*ZF
        ENDIF

!* Rain-out in Cloud
!* relative water gain due to condensation
        IF (ZCLTOTI > ZEPSQLIQ .AND. ZBETASI > 0._JPRB .AND. YAERO_DESC(JT)%RSCAVIN > 0.0_JPRB) THEN
! ice cloud
          ZKNEW = ZKMINICE + ZBETASI/ZCLTOTI
          ZF=ZCLCOV * ZBETASI/ (ZKNEW * ZCLTOTI)
          ZSCAV=YAERO_DESC(JT)%RSCAVIN*ZKNEW*0.5_JPRB
          ZFUNC=ZFUNC+(EXP(-ZSCAV*PTSPHY)-1._JPRB)*ZF
        ENDIF

        IF (.NOT. LLWDAER .AND. ZCLTOT > ZEPSQLIQ .AND. ZBETA > ZEPSFLX) THEN

          ZBETAR=RG*ZBETA / (ZDP * ZCLTOT) ! in s-1
!* safety check
          ZBETAR=MAX( ZBETAR, 0._JPRB )
!* scavenging coefficient  in s-1 
 !In-cloud scavenging is different for aerosol than for gas-phase 
          ZSCAV=(ZRET*ZLIQ2TOT + ZICE2TOT)*ZBETAR
          ZFUNC=ZFUNC+(EXP(-ZSCAV*PTSPHY)-1._JPRB) * ZCLCOV
        ENDIF

      ZDX = (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY)*ZFUNC      ! in kg kg-1
      PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + (ZDX)/PTSPHY                ! in kgkg-1 s-1
      ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - (ZDX)/PTSPHY * ZDP/RG ! in kg m-2 s-1


      !* reevaporation  ! bottom flux smaller and top flux above limit
      IF (ZBETA < 0.0_JPRB .AND. ZFLXR+ZFLXS > ZEPSFLX) THEN
        ZBETAR=ZBETA/(ZFLXR+ZFLXS)
        IF ( ZFLXRB+ZFLXSB  <= ZEPSFLX) THEN
!-- total reevaporation, bottom flux <= ZEPSFLX
         ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)              ! N.D.
       ELSE
!-- partial reevaporation
         ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)        ! N.D.
        ! ZBETAR=ZBETAR*ZFRAC
         ZBETAR=ZBETAR*(1._JPRB-EXP(-2._JPRB*ZBETAR**0.5_JPRB)* &
             & (1._JPRB + 2._JPRB*ZBETAR**0.5_JPRB + 2._JPRB*ZBETAR +1.33_JPRB*ZBETAR**1.5))* &
             & (1._JPRB-ZBETAR) + ZBETAR**2._JPRB

       ENDIF
       ZDX = ZBETAR*ZWLOSS(JL,IWETDEP) * (RG/ZDP)*PTSPHY

       PTENC1(JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY
       ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG ! in kg m-2 s-1
     ENDIF



    !!*   wash-out with rain  when flux is positive and no rain formation or
    !evaporation
!    IF (ZFLXR > ZEPSFLX .AND. ABS(ZBETA) <= ZEPSFLX ) THEN
!*   wash-out with rain  when flux is positive
    IF ((ZFLXR+ZFLXS) > ZEPSFLX) THEN
!    IF (ZCLTOT > ZEPSQLIQ .AND. (ZFLXR+ZFLXS) > ZEPSFLX) THEN
       IF ( IMODE == 1 ) THEN
! calculate fraction in rainwater according to Henry
!* rain water in box in kg/kg
!          ZRAINW=(ZFLXR*PTSPHY)*RG /ZDP   ! eq 16 .. in kg/kg is now input to
!          routine
           IF (LLCONV ) THEN
             ZRAINW=(ZFLXR/(ZFALLSP*ZRHO)) / ZPRCOV   ! that's is better
           ELSE
             ZRAINW = (PRAIN(JL,JK)+PSNOW(JL,JK)) / ZPRCOV ! + PSNOW(JL,JK) !use prognostic variable from input
           ENDIF

          IF ( LLWDAER ) THEN
            ! Rain:
            !ZSCAV = ZFLXR * 1E-3
            ZSCAV = ZFLXR *  YAERO_DESC(JT)%RSCAVBCR
            ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ZFUNC=MAX(0.75_JPRB,ZFUNC)
            ZDX = (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1

            ! add for snow: assume 5E-3 m2/kg
!            ZSCAV = ZFLXS * 5.E-3_JPRB
            ZSCAV = ZFLXS *  YAERO_DESC(JT)%RSCAVBCS
            ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ZFUNC=MAX(0.75_JPRB,ZFUNC)
            ZDX = ZDX + (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1

          ELSE

! gas phase removal
 !* liquid phase to gas phase fraction - henry's law
          ZLIQ2GAS=ZHNRYEFT*ZRAINW*ZRD*PTP(JL,JK)*ZLMMR2VMR
!* liquid phase to total fraction for Henry equilibrium ! used if Henry limited
          ZLIQ2TOT=ZLIQ2GAS/(1.0_JPRB+ZLIQ2GAS)  ! eq 15.

!* rain out limted by mass transfer (Aerosols and HNO3)
          ZSCAV=ZFLXR*0.1_JPRB ! exponent of eq 14, k'i=1 cm-1 = 0.1 mm-1

          ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
! compare henry limited vs mass transfer limited
            IF ( ZLIQ2TOT > (1.0_JPRB -  ZFUNC)  ) THEN
              ZDX = (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1
            ELSE
!  dissolved fraction is flux from above (ZWLOSS) + disolved rain water content
!  ! eq 17  minus delta m_i_top
              ZDX = MIN(0.0_JPRB, (1.0_JPRB - ZLIQ2TOT) * ( ZWLOSS(JL,IWETDEP) * PTSPHY) / (ZDP/RG) - ZLIQ2TOT *&
                & (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY) * ZPRCOV )   ! in kg kg-1
            ENDIF
          ENDIF
        PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
        ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG ! in kg m-2 s-1

! TM5 mode
         ELSEIF (IMODE == 2 ) THEN
           !Roelofs and Lelieveld, 1995
           ! vol mixing ratio
           ZRLWC  = 72._JPRB*((ZFLXR*3600)**0.88)*1.E-9_JPRB
           ! Droplet radius [cm]
           ZRDRAD = 0.1_JPRB*0.3659_JPRB*(((ZFLXR+ZFLXS)*3600)**0.21_JPRB)
           !*******************************************
           ! ZRU: Terminal velocity  in cm/s
           ! ZNRE: Reynolds number
           ! ZNSH: Sherwood number
           ! ZNSC: Schmidt number
           !*******************************************
           ZRU    = 100.*9.58*(1.-EXP(-(ZRDRAD*10./0.885)**1.147))
           ! see Seinfeld (1986)
           ZNRE  = 2.*ZRDRAD*ZRU/ZDGAIR
!!           ZNSC  = (ZDGAIR/ZDGHNO3)**(1./3.)
           ZNSH  = 1.+0.3*(ZNRE**0.5)*ZNSC
           ZKG   = ZDGHNO3/ZRDRAD*ZNSH
           ! this loss rate corresponds to HNO3
           ZBETAR= 3.*ZKG*ZRLWC/ZRDRAD
           ! Adapt for any species
           ZBETAR = ZBETAR * SQRT(ZXMHNO3/YCHEM(JT)%RMOLMASS  )    
!* scavenging rate 
           ZSCAV= ZLIQ2TOT*ZBETAR 
           IF (LLWDAER ) THEN
! Assume a wash-out coefficient of 0.05 mm^-1 (raindepth)
           ! (Dana and Hales, Atmos. Env. 1976, pp. 45-50)
! Note that units of ZFLXR (kg/m2/s) corresponds to mm / s 
             ZSCAV = 0.05 * ZFLXR
             ZSCAV = MAX(ZSCAV, 0._JPRB)
           ENDIF
           ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
           ZDX = (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1
           PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
           ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1
         ENDIF 
     ENDIF  ! rainout
    ENDDO
  ENDDO
ENDDO
 !* LOOP OVER LAYERS  
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PFAERO(JL,:) = ZWLOSS(JL,:)
  ENDDO
ENDDO


!! debug undo anything  
! PTENC1(KIDIA:KFDIA,1:KLEV,1:KCHEM) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KCHEM)
!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AERC_SCAVL19',1,ZHOOK_HANDLE)
END SUBROUTINE AERC_SCAVL19
