! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CHEM_SCAV &
 & ( YDCHEM,YDECLDP,YDEDBUG,YGFL,KIDIA , KFDIA  , KLON , KLEV , KCHEM , KCHEM_SCAV, KSTEP, KCHEM_WETDEP, PTSPHY, &
 &   PRSF1, PDP , PTP, PFLXR, PFLXS, PCLCOV, PCLWAT, PCLICE, PRAIN, PSNOW, PLSM, PCEN, PTENC0, &
 &   PTENC1, PWLOSS, PPRCOV  )

!*** * CHEM_SCAV* - IN-CLOUD AND BELOW CLOUD SCAVENGING OF TRACERS
!      CALLED SEPARATELY FOR CONVECTIVE AND LARGE-SCALE PRECIP. 
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
! PSNWO (KLON,KLEV)           :  Snow water content    (kg/kg) for stratiform precip
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
! PLOSS(KLON,KCHEM_SCAV)          : Total Mass Loss due to scavening     (kg/m2 s-1)!
!
!**   INTERFACE.
!     ----------
!          *CHEM_SCAV* IS CALLED FROM *CHEM_MAIN*.
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
!        S. Remy  : add option for Luo19 scavenging rate 2021-9-8, different
!        values for H+ over ocean and continent, re-evaporation following De
!        Bruine et al. 2018, GMD
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD, RLVTT ,RLSTT ,RTT 
USE YOMLUN   , ONLY : NULOUT
USE YOEDBUG  , ONLY : TEDBUG
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOECLDP  , ONLY : TECLDP
USE YOETHF   , ONLY :  R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,&
 & RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
USE YOMCHEM  , ONLY : TCHEM
IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

TYPE(TCHEM)       ,INTENT(INOUT):: YDCHEM
TYPE(TECLDP)      ,INTENT(INOUT):: YDECLDP
TYPE(TEDBUG)      ,INTENT(INOUT):: YDEDBUG
TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KCHEM, KCHEM_SCAV, KSTEP, KCHEM_WETDEP

REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV) , PRSF1(KLON,KLEV) , PTP(KLON,KLEV)   
REAL(KIND=JPRB),INTENT(IN)    :: PCLCOV(KLON,KLEV) , PCLWAT(KLON,KLEV),PCLICE(KLON,KLEV) , PRAIN(KLON,KLEV), PSNOW(KLON,KLEV)  
REAL(KIND=JPRB),INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB),OPTIONAL, INTENT(IN)    ::  PPRCOV(KLON,KLEV)  
REAL(KIND=JPRB),INTENT(IN)    :: PFLXR(KLON,KLEV+1), PFLXS(KLON,KLEV+1)
REAL(KIND=JPRB),INTENT(IN)    :: PTENC0(KLON,KLEV,KCHEM), PCEN(KLON,KLEV,KCHEM)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY

REAL(KIND=JPRB),INTENT(OUT)   :: PWLOSS(KLON,KCHEM_SCAV), PTENC1(KLON,KLEV,KCHEM)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM), PARAMETER :: IMODE = 1 ! Choice of parameter scheme. Scheme 2 is optimized for TM5
INTEGER(KIND=JPIM) :: JK, JL, JT, IWETDEP

! indices for various chemical species (used to speed up comparisons)
INTEGER(KIND=JPIM), PARAMETER :: IAERO = 1_JPIM
INTEGER(KIND=JPIM), PARAMETER :: IH2O2 = 2_JPIM
INTEGER(KIND=JPIM), PARAMETER :: IHNO3 = 3_JPIM
INTEGER(KIND=JPIM), PARAMETER :: INH3 = 4_JPIM
INTEGER(KIND=JPIM), PARAMETER :: ISO2 = 5_JPIM

! "default" index value (for species not used in comparisons)
INTEGER(KIND=JPIM), PARAMETER :: INOTUSED = -2_JPIM
INTEGER(KIND=JPIM), DIMENSION(KCHEM) :: ICHEM

REAL(KIND=JPRB) :: ZSCAV, ZSCAV_COEF
REAL(KIND=JPRB) :: ZBETA,  ZBETAR, ZBETARI, ZBETASI
REAL(KIND=JPRB) :: ZDP, ZDX, ZDZ, ZFRAC, ZFUNC, ZRHO
REAL(KIND=JPRB) :: ZEPSQLIQ, ZEPSFLX
REAL(KIND=JPRB) :: ZRD, ZLMMR2VMR, ZH2R, ZRET, ZHNRYEFT, ZHNRYEF    
REAL(KIND=JPRB) :: ZLIQ2TOT, ZICE2TOT, ZLIQ2GAS, ZRAINW, ZFALLSP
REAL(KIND=JPRD) :: ZICE2GAS
REAL(KIND=JPRB) :: ZPRCOV, ZCLCOV, ZFLXR, ZFLXS,ZFLXRB, ZFLXSB,  ZCLWAT, ZCLICE,ZCLTOT, ZMINCLCOV, ZMINPRCOV
REAL(KIND=JPRB) :: ZCLWAT2, ZCLICE2,ZCLTOTW,ZCLTOTI, ZF, ZKNEW, ZKMIN,ZKMINICE, ZCLTOT2

 ! Interstitial Fraction: 30% of aerosol remains in atmosphere
!obsolete REAL(KIND=JPRB)  ::  ZINTERST_FR  
REAL(KIND=JPRB),PARAMETER  :: ZDGHNO3 = 0.136      ! viscosity of HNO3 in [cm2/s] 
REAL(KIND=JPRB),PARAMETER  :: ZDGAIR  = 0.133      ! viscosity of air in [cm2/s] 
REAL(KIND=JPRB),PARAMETER  :: ZXMHNO3    =1.008_JPRB + 14.007_JPRB + 3*16.0_JPRB ! HNO3 tracer mass
REAL(KIND=JPRB)           :: ZHPLUS
REAL(KIND=JPRB),PARAMETER :: ZHPLUS_OCEAN =2.51188E-6_JPRB ! [H+] corresponding to rain water pH=5.6 over oceans
REAL(KIND=JPRB),PARAMETER :: ZHPLUS_LAND =1.E-5_JPRB ! [H+] of rain water pH=5 over land

REAL(KIND=JPRB)  :: ZRLWC,ZRDRAD, ZRU, ZNRE,ZNSC,ZNSH,  ZKG, ZTR, ZKSO2, ZKHSO3, ZFACTSO2
REAL(KIND=JPRB)  :: ZHNH3, ZKA1,ZKW
REAL(KIND=JPRB)            :: ZRL      ! composite factor of Rgas and liquid water content of raining cloud
                                       ! rgas (8.314 J/mol/K) ---> 0.08314 atm/(mol/l)/K
                                       ! 1e-6 corresponds to 1 g/m3 dimensionless 
LOGICAL :: LLWDAER 
LOGICAL :: LLPRINT, LLCHEM_WDFR, LLCONV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "abor1.intfb.h"
!#include "fccld.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_SCAV',0,ZHOOK_HANDLE)
ASSOCIATE(YCHEM=>YGFL%YCHEM, &
 & LCHEM_WDFR=>YDCHEM%LCHEM_WDFR, &
 & RCOVPMIN=>YDECLDP%RCOVPMIN, &
 & KSTPDBG=>YDEDBUG%KSTPDBG, NSTPDBG=>YDEDBUG%NSTPDBG)
LLPRINT=.FALSE.
LLCHEM_WDFR=LCHEM_WDFR
LLCONV=.FALSE.

!VH ZINTERST_FR = 0.3_JPRB
!VH aerosol wet deposition in this routine appears too efficient.
!VH instead use general scavening coefficient SCAV_COEF according to Bey (JGR 2011)
  
! deposition for convective Precip
IF (.NOT. PRESENT(PPRCOV) )  LLCONV=.TRUE.

! init indices for chemical species in the YCHEM array (to be used later)
DO JT=1,KCHEM
  SELECT CASE (TRIM(YCHEM(JT)%CNAME))
    CASE ("H2O2")
      ICHEM(JT) = IH2O2
    CASE ("HNO3")
      ICHEM(JT) = IHNO3
    CASE ("NH3")
      ICHEM(JT) = INH3
    CASE ("SO2")
      ICHEM(JT) = ISO2
    CASE ("SO4", "NH4", "Pb", "NO3_A", "PM10", "PM25")
      ICHEM(JT) = IAERO
    CASE DEFAULT
      ICHEM(JT) = INOTUSED 
  END SELECT
ENDDO

!* set flux to zero 
PWLOSS(:,:)=0._JPRB

DO JL=1,NSTPDBG
  IF (KSTEP == KSTPDBG(JL)) THEN
    LLPRINT=.TRUE.
  ENDIF
ENDDO

ZEPSFLX =1.E-18_JPRB
ZEPSQLIQ=1.E-18_JPRB

!* mini cloud cover and precip cover
ZMINCLCOV=0.001_JPRB
ZMINPRCOV=0.001_JPRB

! For L19 rainout
ZKMIN=1.E-4_JPRB
ZKMINICE=1.E-6_JPRB

!* precip fall speed 
ZFALLSP = 5.0_JPRB

!* ideal gas constant in atm M-1
ZRD=1000.0_JPRB * RD * 9.8692_JPRB / 1000000.0_JPRB
ZRD=0.082_JPRB
ZRL=RG/1E2_JPRB*1E-6_JPRB 

!* initialisation    
!* re-evaporation fraction to account for drop shrinking without releasing species
! Jacob says 
!ZFRAC=0.5_JPRB
ZFRAC=0.2_JPRB

!* update tendecies  

!* if not LLCHEM_WDFR 
ZPRCOV =  1.0_JPRB
ZCLCOV =  1.0_JPRB

PTENC1(KIDIA:KFDIA,1:KLEV,1:KCHEM) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KCHEM)

!* LOOP OVER LAYERS  
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (PLSM(JL) >= 0.1_JPRB) THEN
      ZHPLUS=ZHPLUS_LAND
    ELSE
      ZHPLUS=ZHPLUS_OCEAN
    ENDIF
!- precip flux at box top 
    ZFLXR =  PFLXR(JL,JK)
    ZFLXS =  PFLXS(JL,JK)
    ZFLXRB =  PFLXR(JL,JK+1)
    ZFLXSB =  PFLXS(JL,JK+1)

!--  Precip flux change grid box average - indicator for rain/snow formation or evaporation
    ZBETA = ( ZFLXRB - ZFLXR ) + ( ZFLXSB - ZFLXS  )  ! in kg m-2 s-1 
    ZBETARI = ZFLXRB - ZFLXR
    ZBETASI = ZFLXSB - ZFLXS
    ZCLWAT =  PCLWAT(JL,JK)
    ZCLICE =  PCLICE(JL,JK)
    ZCLWAT2 =  MAX(0.0_JPRB,PCLWAT(JL,JK))
    ZCLICE2 =  MAX(0.0_JPRB,PCLICE(JL,JK))


    IF (  ZBETA ==  0.0_JPRB .AND. ZFLXR < ZEPSFLX  )  CYCLE 
    IF ( LLCHEM_WDFR ) THEN 
!* use effective flux and cloud cover for liquid water content
      IF (.NOT. LLCONV ) THEN
        ZPRCOV = MIN(MAX(RCOVPMIN,PPRCOV(JL,JK)),1.0_JPRB)  ! grid scale precip
      ELSE
        ZPRCOV = 0.05    ! convective precip  
        ZPRCOV = 0.1    ! convective precip  
      ENDIF  
      ZCLCOV = MIN(1.0_JPRB,MAX(0.0_JPRB,PCLCOV(JL,JK)))
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
      IF (ZCLCOV >= ZMINCLCOV ) THEN 
        ZCLWAT =  MAX(0.0_JPRB,PCLWAT(JL,JK) /  ZCLCOV)
        ZCLICE =  MAX(0.0_JPRB,PCLICE(JL,JK) /  ZCLCOV)
      ELSE
        ZCLWAT = 0.0_JPRB  
        ZCLICE = 0.0_JPRB  
      ENDIF 
    ENDIF

    ZCLTOT=ZCLWAT + ZCLICE 
!* calculate air density, layer depth 
    ZRHO=PRSF1(JL,JK)/(RD*PTP(JL,JK))
    ZDZ= PDP(JL,JK) / (ZRHO*RG)
    ZDP =  PDP(JL,JK)
    ZBETARI=ZBETARI*RG/ZDP  ! in s-1
    ZBETASI=ZBETASI*RG/ZDP  ! in s-1
!* total cloud and ice water content
    ZCLTOTW=ZCLWAT2 + PTSPHY * ZBETARI
    ZCLTOTI=ZCLICE2 + PTSPHY * ZBETASI
    ZCLTOT2=ZCLTOTW+ZCLTOTI

! clwc mass mixing ratio to volume mixing ratio : rho air / rho cloud water - note it is liquid water not vapor
    ZLMMR2VMR=ZRHO/1000.0_JPRB 
    IWETDEP=0

!* LOOP over species
    DO JT=1,KCHEM

       IF ( YCHEM(JT)%HENRYA <=  0.0_JPRB ) CYCLE

! identify aerosols
       LLWDAER = (ICHEM(JT) == IAERO) 

       IWETDEP=IWETDEP+1
       ZHNRYEF=YCHEM(JT)%HENRYA
       ZH2R  =YCHEM(JT)%HENRYB

       IF ( LLWDAER ) THEN 
! Introduce scavenging coefficients for aerosol according to
! Bourgeois and Bey, JGR 2011
         IF (PTP(JL,JK) < 248.0_JPRB) THEN
            ! for stratiform ice clouds:
            ZSCAV_COEF = 0.06
         ELSEIF (PTP(JL,JK) < 268.0_JPRB) THEN
            ! For stratiform mixed clouds: 
            ZSCAV_COEF = 0.06
            ! Adaptation needed for Coarse Soluable (0.75) and Coarse Insoluable (0.4)
         ELSE ! PTP > 273 Kelvin
            ! For stratiform liquid clouds:
            ZSCAV_COEF = 0.1
            ! Variation needed for different kind of particles, 
            ! ranging from 0.06 (nucl mode)-0.99 ( coarse soluable )
            ! Value of 0.1 is based on tuning - not yet best.
         ENDIF

      ELSE 
! Gas-phase removal, based on Henry solubility


!* scavenging coefficient according to "Harvard wet deposition scheme for GMI" Jacobs, D.J. 2000
! zscav should be solved fraction ->  henry 
!* Retention coefficient for mixed clouds
        IF (PTP(JL,JK) < 268.0_JPRB) THEN
          IF (ICHEM(JT) == IH2O2) THEN
            !ZRET=0.05_JPRB   ! Special case for H2O2 Jacobs Table 1
            ZRET=0.6_JPRB   ! Use optimized retention, as H2O2 wet dep. appears too low 
                             ! compared to TM5. See also von Blohn et al., ACP 2011
          ELSEIF (ICHEM(JT) == IHNO3 .OR. ICHEM(JT) == INH3) THEN
            ZRET=1.0_JPRB   ! Keep unity (Mari et al.)
          ELSE
            ZRET=0.02_JPRB   ! For others assume retention of 0.02
          ENDIF
        ELSE
!* otherwise 1.0  
          ZRET=1.0_JPRB
        ENDIF

!* Henry coefficicnet 
        ! Henrys are in M/atm = (mol/liter)/atm - RD has to be changed 
        ! Temperature dependency of henry constant 
        ZTR=(1.0_JPRB/PTP(JL,JK)-1.0_JPRB/298.0_JPRB)
        ZHNRYEFT=ZHNRYEF*EXP(ZH2R * ZTR)  
! effective Henry coeff for SO2 (Seinfeld Pandis, 1998. p 350 / Table 6.4)
        IF (ICHEM(JT) == ISO2) THEN
           ZKSO2=1.3E-2_JPRB*EXP(2090.0_JPRB*ZTR)     ! Aq. equi. const. so2<=>hso3m+hplus
           ZKHSO3=6.6E-8_JPRB*EXP(1510.0_JPRB*ZTR)    ! Aq. equi. const. hso3m<=>so3-- + hplus
           ZFACTSO2=1.0 + ZKSO2 / ZHPLUS  + (ZKSO2*ZKHSO3)/(ZHPLUS**2.0_JPRB)
           ZHNRYEFT=ZFACTSO2 * 1.23_JPRB * EXP ( 3000.0_JPRB * ZTR)
        ELSEIF (ICHEM(JT) == INH3) THEN
! Effective Henry coeff for NH3 (Seinfeld and Pandis, p. 353 / Table 6.4)
           ZHNH3=ZHNRYEF*EXP( ZH2R*ZTR )
           ZKA1=1.7E-5_JPRB *EXP( -450._JPRB*ZTR )
           ZKW =1.0E-14_JPRB*EXP( -6718._JPRB*ZTR )
           ZHNRYEFT = ZHNH3*(1._JPRB + ZKA1*ZHPLUS/ZKW)
        ENDIF                   

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
        ZICE2GAS=0.0_JPRD 
        IF (ICHEM(JT) == IHNO3 .AND. ZCLICE > 0.0_JPRB) THEN
          ZICE2GAS=ZHNRYEFT*ZCLICE*ZRD*PTP(JL,JK)*ZLMMR2VMR
        ENDIF
        IF (ICHEM(JT) == IH2O2 .AND. ZCLICE > 0.0_JPRB) THEN
! calculate following Mari et al. 2000 as decribed in Jacob, eq 9 
!         ZVPICE=FOEEICE(PTP(JL,JK)) ! water vapour saturated over ice in Pa  
!         ZICE2GAS  = ((ZCLICE*ZLMMR2VMR)/(ZVPICE/PRSF1(JL,JK))) *0.82_JPRB ! this seems far to low ????
! calculate following Lawrence and Crutzen 1998, Tellus, as cited in Neu and Prather et al. 2011 (ACP)  
           ZICE2GAS = (ZCLICE*ZLMMR2VMR) * 5.0E4_JPRD * EXP( 0.48_JPRD * & 
           &          10.0_JPRD**( - (PTP(JL,JK) - 273.15_JPRD) / 43.0_JPRD))
!          IF (ZCLICE > 1.0e-4 .AND. FLOOR(PTP(JL,JK)) == 258 ) THEN
!            print*,TRIM(YCHEM(JT)%CNAME), JL, JK, ZCLICE , PTP(JL,JK) ,  ZICE2GAS , ZTEST, ZTEST2, ZVPICE, PRSF1(JL,JK)
!          ENDIF  
       ENDIF
!* liquid phase to total fraction 
        ZLIQ2TOT=ZLIQ2GAS/(1.0_JPRB+ZLIQ2GAS+ZICE2GAS)  
!* ice phase to total fraction 
        ZICE2TOT=ZICE2GAS/(1.0_JPRB+ZLIQ2GAS+ZICE2GAS)  

      ENDIF ! Definition of henry-solubility (LWDAER=FALSE) 

      SELECT CASE (KCHEM_WETDEP )

        CASE (1) ! original parameterization

!* Rain-out in Cloud
          IF (ZCLTOT > ZEPSQLIQ .AND. ZBETA > 0.0_JPRB ) THEN
!* relative water loss due to precip formation
            ZBETAR=ZBETA / ( ZDP * ZCLTOT / RG ) ! in s-1
            IF (LLPRINT .AND. ZBETA /= 0._JPRB .AND. ZCLWAT > ZEPSQLIQ ) THEN
              WRITE(UNIT=NULOUT,FMT='(1x,''ZSCAV1'',5I3,(8E10.3))') KSTEP,JL,JK,&
&             ZBETAR,ZBETA,ZSCAV,ZDZ,ZRHO,PCLCOV(JL,JK),ZCLWAT
            ENDIF

!* safety check
            ZBETAR=MAX( ZBETAR, 0._JPRB )

!* scavenging coefficient  in s-1 
 !In-cloud scavenging is different for aerosol than for gas-phase 
            IF ( LLWDAER ) THEN
         ! aerosol: no longer use interstitial fraction 
          ! ZSCAV=(1.0_JPRB-ZINTERST_FR)*ZBETAR
               ZSCAV = ZSCAV_COEF * ZBETAR
            ELSE
              ZSCAV=(ZRET*ZLIQ2TOT + ZICE2TOT)*ZBETAR
            ENDIF  
!
            ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ZDX = PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB) * ZCLCOV      ! in kg kg-1
            PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
            PWLOSS(JL,IWETDEP) = PWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1
            IF (LLPRINT .AND. ZBETA > 0._JPRB) THEN
              WRITE(UNIT=NULOUT,FMT='(1x,''ZSCAV2'',5I3,(9E10.3))')&
&             KSTEP,JL,JK,PCEN(JL,JK,JT),PTENC1(JL,JK,JT),PTENC0(JL,JK,JT),ZBETA,ZFUNC,ZDX,PWLOSS(JL,IWETDEP)
            ENDIF
          ENDIF

        CASE (2) ! Luo et al. 2019 scavenging rates

          ZSCAV=0._JPRB
          ZFUNC=0._JPRB
!* Rain-out in Cloud
          IF (ZCLTOTW > ZEPSQLIQ .AND. ZBETARI > 0._JPRB) THEN
!* relative water loss due to precip formation

        !* relative water loss due to precipA  formation
            ZKNEW = ZKMIN + ZBETARI/ZCLTOTW
            ZF=ZCLCOV * ZBETARI/ (ZKNEW * ZCLTOTW)
            IF ( LLWDAER ) THEN
              ZSCAV=ZKNEW
            ELSE
              ZSCAV=ZKNEW*ZRET*ZLIQ2TOT
            ENDIF
            ZFUNC=(EXP(-ZSCAV*PTSPHY)-1._JPRB)*ZF
            IF (LLPRINT) THEN
              WRITE(*,*) "ZSCAV1_LIQ", YCHEM(JT)%CNAME,JK,&
& ZBETARI,ZKNEW,ZF,ZSCAV,ZFUNC,ZDZ,ZRHO,PCLCOV(JL,JK),ZCLWAT,ZLIQ2TOT,ZLIQ2GAS,ZICE2GAS,ZCLTOTW
            ENDIF
          ENDIF

!* Rain-out in Cloud
!* relative water gain due to condensation
          IF (ZCLTOTI > ZEPSQLIQ .AND. ZBETASI > 0._JPRB) THEN
! ice cloud
            ZKNEW = ZKMINICE + ZBETASI/ZCLTOTI
            ZF=ZCLCOV * ZBETASI/ (ZKNEW * ZCLTOTI)
            IF ( LLWDAER ) THEN
              ZSCAV=ZKNEW*0.5_JPRB
             ELSE
               ZSCAV=ZKNEW*ZICE2TOT
             ENDIF
             ZFUNC= ZFUNC+(EXP(-ZSCAV*PTSPHY)-1._JPRB)*ZF
             IF (LLPRINT) THEN
               WRITE(*,*) "ZSCAV1_ICE", YCHEM(JT)%CNAME,KSTEP,JL,JK,&
& ZBETASI,ZKNEW,ZF,ZSCAV,ZFUNC,ZDZ,ZRHO,PCLCOV(JL,JK),ZCLWAT,ZICE2TOT,ZLIQ2GAS,ZICE2GAS,ZCLTOTI
            ENDIF
          ENDIF
          IF (ZFUNC < 0._JPRB) THEN
            ZDX = (PCEN(JL,JK,JT)+PTENC0(JL,JK,JT)*PTSPHY)*ZFUNC      ! in kg kg-1
            PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
            PWLOSS(JL,IWETDEP) = PWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG ! in kg m-2 s-1
          ENDIF
        CASE DEFAULT
          CALL ABOR1(" KCHEM_WETDEP option doesn't exist")
      END SELECT
!* reevaporation  ! bottom flux smaller and top flux above limit
      IF (ZBETA < 0.0_JPRB .AND. ZFLXR+ZFLXS > ZEPSFLX ) THEN
        ZBETAR=ZBETA/(ZFLXR+ZFLXS)
        IF ( ZFLXRB+ZFLXSB  <= ZEPSFLX) THEN
!-- total reevaporation, bottom flux <= ZEPSFLX
         ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)              ! N.D.
       ELSE
         !-- partial reevaporation, following De Bruine et al. 2018, GMD
         ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)        ! N.D.
         ZBETAR=ZBETAR*(1._JPRB-EXP(-2._JPRB*ZBETAR**0.5_JPRB)* &
             & (1._JPRB + 2._JPRB*ZBETAR**0.5_JPRB + 2._JPRB*ZBETAR+1.33_JPRB*ZBETAR**1.5))* &
             & (1._JPRB-ZBETAR) + ZBETAR**2._JPRB
       ENDIF
       ZDX = (ZBETAR*PWLOSS(JL,IWETDEP) * (RG/ZDP) *PTSPHY) * ZPRCOV 
       PTENC1(JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY
       PWLOSS(JL,IWETDEP) = PWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1
       IF (LLPRINT .AND. ZBETAR > 0._JPRB) THEN
         WRITE(UNIT=NULOUT,FMT='(1x,''ZSCAV3'',5I3,(11E10.3))') KSTEP,JL,JK,&
&        PTENC0(JL,JK,JT),PCEN(JL,JK,JT),PTENC1(JL,JK,JT),PWLOSS(JL,IWETDEP),ZBETAR,ZDX
       ENDIF
     ENDIF

!!*   wash-out with rain  when flux is positive and no rain formation or evaporation 
!    IF (ZFLXR > ZEPSFLX .AND. ABS(ZBETA) <= ZEPSFLX ) THEN    
!*   wash-out with rain  when flux is positive  
    IF (ZFLXR > ZEPSFLX .OR. ZFLXS > ZEPSFLX) THEN    
       IF ( IMODE == 1 ) THEN
! calculate fraction in rainwater according to Henry 
!* rain water in box in kg/kg
!          ZRAINW=(ZFLXR*PTSPHY)*RG /ZDP   ! eq 16 .. in kg/kg is now input to routine 
           IF (LLCONV ) THEN 
             ZRAINW=(ZFLXR/(ZFALLSP*ZRHO)) / ZPRCOV   ! that's is better
           ELSE        
             ZRAINW = (PRAIN(JL,JK)+PSNOW(JL,JK)) / ZPRCOV ! + PSNOW(JL,JK) ! use prognostic variable from input 
           ENDIF

          IF ( LLWDAER ) THEN
! aerosol removal: rain
            ! ZSCAV=(8.4E-6_JPRB*(ZFLXR*3600_JPRB)**0.79_JPRB ) ! from H.Webster and Thomson,2014
            ! ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ! ZDX = PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1
            ! ZFUNC=MAX(0.97_JPRB, ZFUNC)
! aerosol removal: snow
            ! ZSCAV=(8.E-6_JPRB*(ZFLXS*3600_JPRB)**0.305_JPRB ) ! from H.Webster and Thomson,2014
            ! ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ! ZDX = ZDX + PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1
            ! ZFUNC=MAX(0.97_JPRB, ZFUNC)

            !alternatively, follow Croft et al. (ACP 2009)
            ! where ZCOEF_BC is defined such that this matches with
            ! the assumed rain drop diameter of 4 mm
            ! They provide values ranging from 0.1 (coarse soluable)
            ! to 1e-4 (Aitken).
            ! Here assume aerosol on accumulation mode (1E-3 m^2 kg^(-1))

            ! Rain:
            ZSCAV = ZFLXR * 1E-3 
            ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ZDX = PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1

            ! add for snow: assume 5E-3 m2/kg 
            ZSCAV = ZFLXS * 5.E-3_JPRB 
            ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            ZDX = ZDX + PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1

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
              ZDX = PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1
            ELSE
!  dissolved fraction is flux from above (PWLOSS) + disolved rain water content ! eq 17  minus delta m_i_top      
              ZDX = MIN(0.0_JPRB, (1.0_JPRB - ZLIQ2TOT) * ( PWLOSS(JL,IWETDEP) * PTSPHY) / (ZDP/RG) - ZLIQ2TOT *&
                & PCEN(JL,JK,JT) * ZPRCOV )   ! in kg kg-1
            ENDIF
          ENDIF
        PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
        PWLOSS(JL,IWETDEP) = PWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1

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
           ZNSC  = (ZDGAIR/ZDGHNO3)**(1./3.)
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
           ZDX = PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV   ! in kg kg-1
           PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
           PWLOSS(JL,IWETDEP) = PWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1
         ENDIF 
  
     ENDIF  ! rainout
   
    ENDDO
  ENDDO
ENDDO

!! debug undo anything  
! PTENC1(KIDIA:KFDIA,1:KLEV,1:KCHEM) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KCHEM)
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_SCAV',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_SCAV 
