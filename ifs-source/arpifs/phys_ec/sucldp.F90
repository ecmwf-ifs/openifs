! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUCLDP(YDSTA,YDDIMV,YDPHY2,YDECLDP)

!**** *SUCLDP*   - INITIALIZE COMMON YOECLD CONTROLLING *CLOUDSC*

!     PURPOSE.
!     --------
!           INITIALIZE YOECLDP

!**   INTERFACE.
!     ----------
!        CALL *SUCLDP* FROM *SUPHEC*
!              ------        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLDP

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-02-07
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Ko"hler &   03-Dec-2004 total water variance setup for
!        A.Tompkins                moist advection-diffusion PBL and
!                                  7200s autoconversion timescale (RKCONV)
!                                  instead of 10000s
!        R.Forbes      28-May-2008 Changed factor in RTAUMEL from 1.5 to 0.66 
!        A.Tompkins/JJM 20080729   cloud-aerosol interactions
!        P.Bechtold    19-Jan-2009 Changed RCLDIFF from 3.E-6 to 5.E-6
!        A.Tompkins/RF  Jul 2009   Many new constants for multi-phase microphysics
!                                  Addition of namcldp namelist
!        R.Forbes       20110301   Added ice deposition parameters
!        R.Forbes       20111001   Decreased ice fall speed RVICE
!        N.Semane+P.Becht 20120607 Small planet modified fall speeds for gravity
!        N.Semane+P.Bechtold     04-10-2012 Add RPLRG/RPLDARE/RVRFACTOR factor for small planet
!        R.Forbes      01-Mar-2013 Included parameters for autoconv/acc/evap
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        R.Forbes      15-Jan-2015 Included parameters for ice, snow and rain freezing
!        R.Forbes      15-Dec-2015 Included parameters for rain fallspeed
!        R.Forbes      15-Apr-2017 Removed erroneous lines for RCL_LAM1R,RCL_LAM2R
!                                  and modified turbulent erosion coefficients
!        R.Forbes      Jan-2019    Added new microphys params, can be set in namelist
!                                  Added parameters for new accretion formulation
!        F.Vana        14-Sep-2020 Weight factor for cloud to help adjoint accuracy
!        R.Forbes      Nov-2020    Added RSSICEFACTOR, various constants for ice
!                                  Modified RCLDIFF/RCLCRITSNOW/RVICE/RCL_OVERLAPLIQICE
!        R.Forbes      Dec-2020    Added RCL_LAM1S/2S for snow PSD slope
!        R.Forbes      May-2022    Addded precip type severity definition
!     ------------------------------------------------------------------

USE YOMSTA    , ONLY : TSTA
USE YOMDIMV   , ONLY : TDIMV
USE YOMPHY2   , ONLY : TPHY2
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RTT, RPI
USE YOECLDP   , ONLY : TECLDP
USE YOMLUN    , ONLY : NULNAM, NULOUT  
USE YOMDYNCORE, ONLY : RPLRG, RPLDARE

IMPLICIT NONE

TYPE(TSTA)  ,INTENT(IN) :: YDSTA
TYPE(TDIMV) ,INTENT(IN) :: YDDIMV
TYPE(TPHY2) ,INTENT(IN) :: YDPHY2
TYPE(TECLDP),INTENT(INOUT), TARGET :: YDECLDP

REAL(KIND=JPRB), EXTERNAL :: FCGENERALIZED_GAMMA
REAL(KIND=JPRB) :: ZX, ZPLRG
INTEGER(KIND=JPIM) :: IX, JLEV
REAL(KIND=JPRB) :: ZGAMMA1R,ZGAMMA2R,ZGAMMA3R,ZGAMMA4R,ZGAMMA5R,ZGAMMA6R
REAL(KIND=JPRB) :: ZGAMMA1S,ZGAMMA2S,ZGAMMA3S,ZGAMMA4S
REAL(KIND=JPRB) :: ZGAMMA1I,ZGAMMA2I,ZGAMMA3I,ZGAMMA4I

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

INTEGER(KIND=JPIM), POINTER :: NAERCLD, NAECLBC, NAECLDU, NAECLOM, NAECLSS, NAECLSU
REAL(KIND=JPRB), POINTER :: RCCNOM, RCCNSS, RCCNSU, RVICE, RCLDIFF, RCLDIFF_CONVI, RKCONV
REAL(KIND=JPRB), POINTER :: RKOOPTAU, RLCRITSNOW, RTAUMEL, RDEPLIQREFDEPTH, RDEPLIQREFRATE
REAL(KIND=JPRB), POINTER :: RCL_INHOMOGAUT, RCL_INHOMOGACC, RCL_OVERLAPLIQICE, RCL_EFFRIME
REAL(KIND=JPRB), POINTER :: RTAU_CLD_TLAD
LOGICAL, POINTER :: LCLOUD_INHOMOG
LOGICAL, POINTER :: LCLDBUDC, LCLDBUDL, LCLDBUDI, LCLDBUDT
LOGICAL, POINTER :: LCLDBUD_VERTINT, LCLDBUD_TIMEINT

#include "namcldp.nam.h"

!------------------------
!*       1.    SET VALUES
!------------------------

IF (LHOOK) CALL DR_HOOK('SUCLDP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, TSPHY=>YDPHY2%TSPHY, &
 & LAERICEAUTO=>YDECLDP%LAERICEAUTO, LAERICESED=>YDECLDP%LAERICESED, &
 & LAERLIQAUTOCP=>YDECLDP%LAERLIQAUTOCP, LAERLIQAUTOCPB=>YDECLDP%LAERLIQAUTOCPB, &
 & LAERLIQAUTOLSP=>YDECLDP%LAERLIQAUTOLSP, LAERLIQCOLL=>YDECLDP%LAERLIQCOLL, &
 & LCLDBUDGET=>YDECLDP%LCLDBUDGET, LCLDEXTRA=>YDECLDP%LCLDEXTRA, &
 & NBETA=>YDECLDP%NBETA, NCLDDIAG=>YDECLDP%NCLDDIAG, NCLDTOP=>YDECLDP%NCLDTOP, &
 & NSHAPEP=>YDECLDP%NSHAPEP, NSHAPEQ=>YDECLDP%NSHAPEQ, NSSOPT=>YDECLDP%NSSOPT, &
 & NPTYPE_SEV2WMO=>YDECLDP%NPTYPE_SEV2WMO, &
 & RAMID=>YDECLDP%RAMID, RAMIN=>YDECLDP%RAMIN, &
 & RCCN=>YDECLDP%RCCN, RCLCRIT=>YDECLDP%RCLCRIT, &
 & RCLCRIT_LAND=>YDECLDP%RCLCRIT_LAND, RCLCRIT_SEA=>YDECLDP%RCLCRIT_SEA, &
 & RCLDMAX=>YDECLDP%RCLDMAX, RCLDTOPCF=>YDECLDP%RCLDTOPCF, &
 & RCLDTOPP=>YDECLDP%RCLDTOPP, RCL_AI=>YDECLDP%RCL_AI, &
 & RCL_APB1=>YDECLDP%RCL_APB1, RCL_APB2=>YDECLDP%RCL_APB2, &
 & RCL_APB3=>YDECLDP%RCL_APB3, RCL_AR=>YDECLDP%RCL_AR, RCL_AS=>YDECLDP%RCL_AS, &
 & RCL_BI=>YDECLDP%RCL_BI, RCL_BR=>YDECLDP%RCL_BR, RCL_BS=>YDECLDP%RCL_BS, &
 & RCL_CDENOM1=>YDECLDP%RCL_CDENOM1, RCL_CDENOM2=>YDECLDP%RCL_CDENOM2, &
 & RCL_CDENOM3=>YDECLDP%RCL_CDENOM3, RCL_CI=>YDECLDP%RCL_CI, &
 & RCL_CONST1I=>YDECLDP%RCL_CONST1I, RCL_CONST1R=>YDECLDP%RCL_CONST1R, &
 & RCL_CONST1S=>YDECLDP%RCL_CONST1S, RCL_CONST2I=>YDECLDP%RCL_CONST2I, &
 & RCL_CONST2R=>YDECLDP%RCL_CONST2R, RCL_CONST2S=>YDECLDP%RCL_CONST2S, &
 & RCL_CONST3I=>YDECLDP%RCL_CONST3I, RCL_CONST3R=>YDECLDP%RCL_CONST3R, &
 & RCL_CONST3S=>YDECLDP%RCL_CONST3S, RCL_CONST4I=>YDECLDP%RCL_CONST4I, &
 & RCL_CONST4R=>YDECLDP%RCL_CONST4R, RCL_CONST4S=>YDECLDP%RCL_CONST4S, &
 & RCL_CONST5I=>YDECLDP%RCL_CONST5I, RCL_CONST5R=>YDECLDP%RCL_CONST5R, &
 & RCL_CONST5S=>YDECLDP%RCL_CONST5S, RCL_CONST6I=>YDECLDP%RCL_CONST6I, &
 & RCL_CONST6R=>YDECLDP%RCL_CONST6R, RCL_CONST6S=>YDECLDP%RCL_CONST6S, &
 & RCL_CONST7R=>YDECLDP%RCL_CONST7R, RCL_CONST7S=>YDECLDP%RCL_CONST7S, & 
 & RCL_CONST8R=>YDECLDP%RCL_CONST8R, RCL_CONST8S=>YDECLDP%RCL_CONST8S, &
 & RCL_CONST9R=>YDECLDP%RCL_CONST9R, RCL_CONST10R=>YDECLDP%RCL_CONST10R, &
 & RCL_CONST7I=>YDECLDP%RCL_CONST7I, &
 & RCL_EFF_RACW=>YDECLDP%RCL_EFF_RACW, &
 & RCL_CR=>YDECLDP%RCL_CR, RCL_CS=>YDECLDP%RCL_CS, RCL_DI=>YDECLDP%RCL_DI, &
 & RCL_DR=>YDECLDP%RCL_DR, RCL_DS=>YDECLDP%RCL_DS, &
 & RCL_DYNVISC=>YDECLDP%RCL_DYNVISC, & 
 & RCL_LAMBDA1I=>YDECLDP%RCL_LAMBDA1I, RCL_LAMBDA2I=>YDECLDP%RCL_LAMBDA2I, &
 & RCL_LAM1R=>YDECLDP%RCL_LAM1R, RCL_LAM2R=>YDECLDP%RCL_LAM2R, &
 & RCL_LAM1R_MP=>YDECLDP%RCL_LAM1R_MP, RCL_LAM2R_MP=>YDECLDP%RCL_LAM2R_MP, &
 & RCL_LAM1S=>YDECLDP%RCL_LAM1S, RCL_LAM2S=>YDECLDP%RCL_LAM2S, &
 & RCL_FZRAB=>YDECLDP%RCL_FZRAB, &
 & RCL_FZRBB=>YDECLDP%RCL_FZRBB, RCL_KA273=>YDECLDP%RCL_KA273, &
 & RCL_KKAAC=>YDECLDP%RCL_KKAAC, RCL_KKAAU=>YDECLDP%RCL_KKAAU, &
 & RCL_KKBAC=>YDECLDP%RCL_KKBAC, RCL_KKBAUN=>YDECLDP%RCL_KKBAUN, &
 & RCL_KKBAUQ=>YDECLDP%RCL_KKBAUQ, &
 & RCL_KK_CLOUD_NUM_LAND=>YDECLDP%RCL_KK_CLOUD_NUM_LAND, &
 & RCL_KK_CLOUD_NUM_SEA=>YDECLDP%RCL_KK_CLOUD_NUM_SEA, &
 & RCL_SCHMIDT=>YDECLDP%RCL_SCHMIDT, RCL_X1I=>YDECLDP%RCL_X1I, &
 & RCL_X1R_MP=>YDECLDP%RCL_X1R_MP, RCL_X2R_MP=>YDECLDP%RCL_X2R_MP, &
 & RCL_X4R_MP=>YDECLDP%RCL_X4R_MP, &
 & RCL_X1R=>YDECLDP%RCL_X1R, RCL_X1S=>YDECLDP%RCL_X1S, RCL_X2I=>YDECLDP%RCL_X2I, &
 & RCL_X2R=>YDECLDP%RCL_X2R, RCL_X2S=>YDECLDP%RCL_X2S, RCL_X3I=>YDECLDP%RCL_X3I, &
 & RCL_X3S=>YDECLDP%RCL_X3S, RCL_X4I=>YDECLDP%RCL_X4I, RCL_X4R=>YDECLDP%RCL_X4R, &
 & RCL_X4S=>YDECLDP%RCL_X4S, RCOVPMIN=>YDECLDP%RCOVPMIN, &
 & RDENSREF=>YDECLDP%RDENSREF, RDENSWAT=>YDECLDP%RDENSWAT, &
 & RICEHI1=>YDECLDP%RICEHI1, &
 & RICEHI2=>YDECLDP%RICEHI2, RICEINIT=>YDECLDP%RICEINIT, &
 & RLMIN=>YDECLDP%RLMIN, RNICE=>YDECLDP%RNICE, RPECONS=>YDECLDP%RPECONS, &
 & RPRC1=>YDECLDP%RPRC1, RPRC2=>YDECLDP%RPRC2, RPRECRHMAX=>YDECLDP%RPRECRHMAX, &
 & RSNOWLIN1=>YDECLDP%RSNOWLIN1, RSNOWLIN2=>YDECLDP%RSNOWLIN2, &
 & RTHOMO=>YDECLDP%RTHOMO, RSSICEFACTOR=>YDECLDP%RSSICEFACTOR, &
 & RVRAIN=>YDECLDP%RVRAIN, RVRFACTOR=>YDECLDP%RVRFACTOR, RVSNOW=>YDECLDP%RVSNOW, &
 & STPRE=>YDSTA%STPRE)
! Autoconversion/accretion (KK 2000)
! & RCL_KKAac, RCL_KKBac, RCL_KKAau, RCL_KKBauq, RCL_KKBaun, &
! & RCL_KK_cloud_num_sea, RCL_KK_cloud_num_land, &
! Ice
!& RCL_AI, RCL_BI, RCL_CI, RCL_DI, RCL_X1I, RCL_X2I, RCL_X3I, RCL_X4I, &
!& RCL_CONST1I, RCL_CONST2I, RCL_CONST3I, RCL_CONST4I, RCL_CONST5I, RCL_CONST6I,&
!& RCL_APB1, RCL_APB2, RCL_APB3, &
! Snow
!& RCL_AS, RCL_BS, RCL_CS, RCL_DS, RCL_X1S, RCL_X2S, RCL_X3S, RCL_X4S, &
!& RCL_CONST1S, RCL_CONST2S, RCL_CONST3S, RCL_CONST4S, RCL_CONST5S, &
!& RCL_CONST6S, RCL_CONST7S, RCL_CONST8S, &
! Rain
!& RDENSWAT, RDENSREF, RCL_AR, RCL_BR, RCL_CR, RCL_DR, &
!& RCL_X1R, RCL_X2R, RCL_X4R, RCL_KA273, &
!& RCL_CDENOM1, RCL_CDENOM2, RCL_CDENOM3, RCL_SCHMIDT, RCL_DYNVISC, &
!& RCL_CONST1R, RCL_CONST2R, RCL_CONST3R, RCL_CONST4R, RCL_CONST5R, RCL_CONST6R, &
!& RCL_LAM1R, RCL_LAM2R, RCL_FZRAB, RCL_FZRBB, &


! Associate pointers for variables in namelist
NAERCLD       => YDECLDP%NAERCLD
NAECLBC       => YDECLDP%NAECLBC
NAECLDU       => YDECLDP%NAECLDU
NAECLOM       => YDECLDP%NAECLOM
NAECLSS       => YDECLDP%NAECLSS
NAECLSU       => YDECLDP%NAECLSU
RCCNOM        => YDECLDP%RCCNOM
RCCNSS        => YDECLDP%RCCNSS
RCCNSU        => YDECLDP%RCCNSU
RVICE         => YDECLDP%RVICE
RCLDIFF       => YDECLDP%RCLDIFF
RCLDIFF_CONVI => YDECLDP%RCLDIFF_CONVI
RKCONV        => YDECLDP%RKCONV
RTAU_CLD_TLAD => YDECLDP%RTAU_CLD_TLAD
LCLOUD_INHOMOG=>YDECLDP%LCLOUD_INHOMOG
RCL_INHOMOGAUT  => YDECLDP%RCL_INHOMOGAUT 
RCL_INHOMOGACC  => YDECLDP%RCL_INHOMOGACC 
RCL_OVERLAPLIQICE => YDECLDP%RCL_OVERLAPLIQICE
RCL_EFFRIME     => YDECLDP%RCL_EFFRIME
RKOOPTAU      => YDECLDP%RKOOPTAU
RLCRITSNOW    => YDECLDP%RLCRITSNOW
RTAUMEL       => YDECLDP%RTAUMEL
RDEPLIQREFDEPTH => YDECLDP%RDEPLIQREFDEPTH
RDEPLIQREFRATE  => YDECLDP%RDEPLIQREFRATE
LCLDBUDC      => YDECLDP%LCLDBUDC
LCLDBUDL      => YDECLDP%LCLDBUDL
LCLDBUDI      => YDECLDP%LCLDBUDI
LCLDBUDT      => YDECLDP%LCLDBUDT
LCLDBUD_VERTINT => YDECLDP%LCLDBUD_VERTINT
LCLDBUD_TIMEINT => YDECLDP%LCLDBUD_TIMEINT

!----------------
! Misc constants 
!----------------
RAMID=0.8_JPRB
RPRC1=100._JPRB
RPRC2=0.5_JPRB
RTAU_CLD_TLAD=0.0_JPRB

!------------------
! safety thresholds
!------------------
RCLDMAX=5.E-3_JPRB
RAMIN=1.E-8_JPRB
RLMIN=1.E-8_JPRB

!---------------------------------------------------------
! Coefficients for cloud edge turbulent erosion
!---------------------------------------------------------
RCLDIFF       = 6.0E-6_JPRB*RPLRG*RPLDARE
RCLDIFF_CONVI = 10.0_JPRB ! increased RCLDIFF in convection by factor RCLDIFF_CONVI

!--------------------------------
! Rain auto conversion constants
!--------------------------------
RCLCRIT      = 4.0E-4_JPRB ! critical autoconversion threshold
RCLCRIT_SEA  = 2.5E-4_JPRB ! critical autoconversion threshold over ocean
RCLCRIT_LAND = 5.5E-4_JPRB ! critical autoconversion threshold over land
RKCONV=RPLRG*RPLDARE/6000._JPRB  ! 1/autoconversion time scale (s)
RCL_INHOMOGAUT = 1.5_JPRB ! Fixed inhomogeneity factor for autoconversion
RCL_INHOMOGACC = 3.0_JPRB ! Fixed inhomogeneity factor for accretion

!--------------------------------
! Snow auto conversion constants
!--------------------------------
! Changing RSNOWLIN1, RSNOWLIN2 makes a bigger difference at high levels 
! where IWC values are low,
RSNOWLIN1=RPLRG*RPLDARE/1000._JPRB ! 1000s = Lin et al. 83
RSNOWLIN2=0.030_JPRB  ! 0.025 = Lin et al. 83
RLCRITSNOW=2.0E-5_JPRB ! critical autoconversion threshold

!-----------------------------------
! Base hydrometeor fall speeds m/s
!-----------------------------------
ZPLRG=RPLDARE*SQRT(RPLRG)
RVICE =0.13_JPRB*ZPLRG
RVRAIN=4.0_JPRB*ZPLRG
RVSNOW=1.0_JPRB*ZPLRG

!---------------------------------------------------
! Used in ice fall-speed modification by air density
! Heymsfield and Iaquinta JAS 2000
!--------------------------------------------------- 
RICEHI1=1.0_JPRB/30000._JPRB
RICEHI2=1.0_JPRB/233._JPRB

!---------------------------------------------------
! Default ice number and liq ccn 
!---------------------------------------------------
RNICE=0.027_JPRB 
RCCN=125.0_JPRB

!---------------------------------------------------
! Default ice deposition variables
!---------------------------------------------------
RCLDTOPCF       = 0.01_JPRB   ! Cloud fraction threshold that defines cloud top 
RDEPLIQREFRATE  = 0.5_JPRB    ! Fraction of deposition rate in cloud top layer
RDEPLIQREFDEPTH = 500.0_JPRB/RPLRG  ! Depth of supercooled liquid water layer (m)
RCL_OVERLAPLIQICE = 0.65_JPRB    ! Overlap assumption for liquid and ice for deposition

!---------------------------------------------------
! Default snow riming
!---------------------------------------------------
RCL_EFFRIME      = 1.0_JPRB   ! Efficiency factor for riming (<1) 

!------------------------------
! initial mass of ice particle
!------------------------------
RICEINIT=1.E-12_JPRB

!---------------------------
! ice microphysics constants
!---------------------------
RTHOMO=RTT-38.0_JPRB ! threshold for homogeneous freezing

!--------------------------
! precipitation evaporation
!--------------------------
RPRECRHMAX=0.7_JPRB  ! Max threshold RH for evaporation for
                     ! a precip coverage of zero 
RPECONS=5.44E-4_JPRB*RPLRG*RPLDARE/RG !evaporation rate coefficient
RVRFACTOR=5.09E-3_JPRB*ZPLRG ! Kessler factor for evaporation (Kessler, 1969): Clear-sky precipitation flux (R)--> R/RVRFACTOR
RCOVPMIN=0.1_JPRB  

!----------------------
! timescale for melting
!----------------------
RTAUMEL=2.0_JPRB*3.6E3_JPRB/(RPLRG*RPLDARE)

!----------------------------------------------------
! Timescale for supersaturation to cause overcast sky
!----------------------------------------------------
RKOOPTAU=3.0_JPRB*3.6E3_JPRB/(RPLRG*RPLDARE)

!--------------------------------
! Tompkins supersaturation scheme
!--------------------------------
NSSOPT=1 

!------------------------------------
! Set reference air and water density
!------------------------------------
RDENSWAT = 1000.0_JPRB  ! kg/m3
RDENSREF = 1.0_JPRB     ! Reference air density

!---------------------------
! Constants for rain
!---------------------------

! Terminal fall speed
! Particle of diameter D, fallspeed vt = c*D^d
RCL_CR = 386.8_JPRB*ZPLRG
RCL_DR = 0.67_JPRB

! Particle of diameter D, mass m = a*D^b
! Volume of a sphere = 4/3 Pi r**3 = Pi/6 D**3
! m = RDENSWAT * Pi/6 D**3
RCL_AR = RDENSWAT*RPI/6._JPRB
RCL_BR = 3.0_JPRB

! Particle size distribution (Abel and Boutle)
RCL_X1R = 0.22_JPRB
RCL_X2R = 2.20_JPRB
RCL_X4R = 0.0_JPRB      ! = mu of gamma distribution 0=gaussian

!Particle size distribution (Marshall and Palmer)
RCL_X1R_MP = 8.0E6_JPRB
RCL_X2R_MP = 0.0_JPRB
RCL_X4R_MP = 0.0_JPRB      ! = mu of gamma distribution 0=gaussian

!Other constants
RCL_LAM1R = RCL_AR*RCL_X1R*FCGENERALIZED_GAMMA(RCL_BR+1._JPRB+RCL_X4R)
RCL_LAM2R = 1.0_JPRB/(RCL_BR+1._JPRB+RCL_X4R-RCL_X2R)

! Equivalent for Marshall-Palmer size distribution
RCL_LAM1R_MP = RCL_AR*RCL_X1R_MP*FCGENERALIZED_GAMMA(RCL_BR+1._JPRB+RCL_X4R_MP)
RCL_LAM2R_MP = 1.0_JPRB/(RCL_BR+1._JPRB+RCL_X4R_MP-RCL_X2R_MP)


!-------------------------------------
! Constants for rain evaporation term
!-------------------------------------
RCL_CDENOM1 = 5.57E11_JPRB  ! =Lv**2./(R*ka_273)
RCL_CDENOM2 = 1.03E8_JPRB   ! =Lv/ka_273
RCL_CDENOM3 = 2.04E2_JPRB   ! =R/(chi*100000)
RCL_KA273   = 2.4E-2_JPRB   ! coeff of conductivity at 0 degC
RCL_SCHMIDT = 0.6_JPRB      ! Schmidt number
RCL_DYNVISC = 1.717E-5_JPRB ! Dynamic viscosity

ZGAMMA1R    = FCGENERALIZED_GAMMA((RCL_DR+5._JPRB)/2._JPRB)
ZGAMMA2R    = FCGENERALIZED_GAMMA(2.0_JPRB + RCL_X4R)
ZGAMMA3R    = FCGENERALIZED_GAMMA(RCL_BR+1.0_JPRB+RCL_X4R)
ZGAMMA4R    = FCGENERALIZED_GAMMA(RCL_DR+RCL_BR+1.0_JPRB+RCL_X4R)
ZGAMMA5R    = FCGENERALIZED_GAMMA(3.0_JPRB + RCL_DR + RCL_X4R)
ZGAMMA6R    = FCGENERALIZED_GAMMA(4.0_JPRB + RCL_BR + RCL_X4R)

RCL_CONST1R = 2._JPRB*RPI*RCL_X1R*RPLRG*RPLDARE
RCL_CONST2R = 0.31_JPRB*RCL_SCHMIDT**(0.3333_JPRB)* &
            & ZGAMMA1R*(RCL_CR/RCL_DYNVISC)**0.5_JPRB
RCL_CONST3R = (RCL_DR+5._JPRB)/2._JPRB - RCL_X2R
RCL_CONST4R = 2._JPRB-RCL_X2R
RCL_CONST7R = RCL_CR*ZGAMMA4R/ZGAMMA3R
RCL_CONST8R = 0.78_JPRB*ZGAMMA2R

!-----------------------------------
! Constants for rain accretion term
!-----------------------------------
RCL_CONST9R  = RPI*0.25_JPRB*RCL_X1R*RCL_CR*ZGAMMA5R
RCL_CONST10R = 3.0_JPRB + RCL_DR - RCL_X2R
RCL_EFF_RACW = 0.7_JPRB ! Efficiency of cloud droplet collection by rain drop

! --------------------------------------------
! Constants for rain autoconversion/accretion (Khairoutdinov and Kogan, 2000)
! --------------------------------------------
RCL_KKAAC  = 67._JPRB   ! s-1
RCL_KKBAC  = 1.15_JPRB
RCL_KKAAU  = 1350._JPRB ! s-1
RCL_KKBAUQ = 2.47_JPRB
RCL_KKBAUN = -1.79_JPRB
RCL_KK_CLOUD_NUM_SEA  =  50._JPRB ! m-3
RCL_KK_CLOUD_NUM_LAND = 300._JPRB ! m-3
! default false: use constant ZEaut, ZEacc
! if true, use regime-dependent FSD to calculate
! ZEaut and ZEacc
LCLOUD_INHOMOG=.FALSE.

!-----------------------------
! Constants for rain freezing
!-----------------------------
RCL_FZRAB = -0.66_JPRB ! negative of Ab constant from Wisneretal(1972) (degC-1)
RCL_FZRBB = 2.E2_JPRB ! Bb constant from Wisneretal(1972) (m-3 s-1)

RCL_CONST5R = RPI/(6._JPRB)*RCL_FZRBB*RCL_AR*RCL_X1R*ZGAMMA6R
RCL_CONST6R = -1._JPRB*(4._JPRB+RCL_BR+RCL_X4R-RCL_X2R)

!---------------------------
! Constants for ice
!---------------------------
! Particle of diameter D, mass m = a*D^b
RCL_AI = 0.069_JPRB ! kg m-2
RCL_BI = 2.0_JPRB

! Particle of diameter D, fallspeed vt = c*D^d
RCL_CI = 16.8_JPRB*ZPLRG
RCL_DI = 0.527_JPRB
!Morrison and Gettelman for cloud ice
RCL_CI = 100. !700.
RCL_DI = 1.0
!RCL_FI = 0.0

! Particle size distribution
RCL_X1I = 2.0E6_JPRB ! m-4
RCL_X2I = 0.0_JPRB
RCL_X3I = 1.0_JPRB
RCL_X4I = 0.0_JPRB

! Gamma functions
ZGAMMA1I = FCGENERALIZED_GAMMA(RCL_BI+1.0_JPRB+RCL_X4I)
ZGAMMA2I = FCGENERALIZED_GAMMA(RCL_X4I+2.0_JPRB)
ZGAMMA3I = FCGENERALIZED_GAMMA((RCL_DI+5.0_JPRB+2.0_JPRB*RCL_X4I)/2.0_JPRB)
ZGAMMA4I = FCGENERALIZED_GAMMA(RCL_DI+RCL_BI+1.0_JPRB+RCL_X4I)

!Other constants
RCL_LAMBDA1I = RCL_AI*RCL_X1I*FCGENERALIZED_GAMMA(RCL_BI+1._JPRB+RCL_X4I)
RCL_LAMBDA2I = 1.0_JPRB/(RCL_BI+1._JPRB+RCL_X4I-RCL_X2I)

!Other constants
RCL_CONST1I = 1.0_JPRB/(RCL_AI*RCL_X1I*ZGAMMA1I) ! lambda numerator
RCL_CONST2I = 1.0_JPRB*RPI*RCL_X1I ! 4*Pi*C where C=0.25*D (Westbrook2008) not D/Pi (WB1999)
RCL_CONST3I = 0.44_JPRB*(0.6_JPRB)**(1.0_JPRB/3.0_JPRB)*RCL_CI**(0.5_JPRB) &
 &             *ZGAMMA3I/(1.717E-5_JPRB)**0.5_JPRB
RCL_CONST4I = (2.0_JPRB+RCL_X4I-RCL_X2I)/(RCL_BI+1.0_JPRB+RCL_X4I-RCL_X2I)
RCL_CONST5I = (5.0_JPRB+RCL_DI+2.0_JPRB*RCL_X4I-2.0_JPRB*RCL_X2I) &
 &             /(2.0_JPRB*(RCL_BI+1._JPRB+RCL_X4I-RCL_X2I)) 
RCL_CONST6I = ZGAMMA2I
RCL_CONST7I = RCL_CI*ZGAMMA4I/ZGAMMA1I

RCL_APB1 = 7.14E11_JPRB
RCL_APB2 = 1.16E8_JPRB
RCL_APB3 = 2.416E2_JPRB

!---------------------------
! Constants for snow
!---------------------------
! Particle of diameter D, mass m = a*D^b
RCL_AS = 0.069_JPRB ! kg m-2
RCL_BS = 2.0_JPRB

! Particle of diameter D, fallspeed vt = c*D^d
RCL_CS = 16.8_JPRB*ZPLRG
RCL_DS = 0.527_JPRB

! Particle size distribution
RCL_X1S = 2.0E6_JPRB ! m-4
RCL_X2S = 0.0_JPRB
RCL_X3S = 1.0_JPRB
RCL_X4S = 0.0_JPRB

! Gamma functions 
ZGAMMA1S = FCGENERALIZED_GAMMA(RCL_BS+1.0_JPRB+RCL_X4S)
ZGAMMA2S = FCGENERALIZED_GAMMA(3.0_JPRB+RCL_DS+RCL_X4S)
ZGAMMA3S = FCGENERALIZED_GAMMA((RCL_DS+5.0+2.0*RCL_X4S)/2.0)!Equiv to ZGAMMA3I
ZGAMMA4S = FCGENERALIZED_GAMMA(RCL_X4S+2.0_JPRB) !Equiv to ZGAMMA2I

! Constants for snow size distribution slope (Lambda)
RCL_LAM1S = RCL_AS*RCL_X1S*FCGENERALIZED_GAMMA(RCL_BS+1.0_JPRB+RCL_X4S)
RCL_LAM2S = 1.0_JPRB/(RCL_BS+1._JPRB+RCL_X4S-RCL_X2S)

! Other constants
RCL_CONST1S = 1.0_JPRB/(RCL_AS*RCL_X1S*ZGAMMA1S)! lambda numerator
RCL_CONST2S = 1.0_JPRB*RPI*RCL_X1S ! 4*Pi*C where C=0.25*D (Westbrook2008) not D/Pi (WB1999)
RCL_CONST3S = 0.44_JPRB*(0.6_JPRB)**(1.0_JPRB/3.0_JPRB)*RCL_CS**(0.5_JPRB) &
 &             *ZGAMMA3S/(1.717E-5_JPRB)**0.5_JPRB
RCL_CONST4S = (2.0_JPRB+RCL_X4S-RCL_X2S)/(RCL_BS+1.0_JPRB+RCL_X4S-RCL_X2S)
RCL_CONST5S = (5.0_JPRB+RCL_DS+2.0_JPRB*RCL_X4S-2.0_JPRB*RCL_X2S) &
 &             /(2.0_JPRB*(RCL_BS+1._JPRB+RCL_X4S-RCL_X2S)) 
RCL_CONST6S = ZGAMMA4S

! Constants for snow riming
RCL_CONST7S = RPI*0.25_JPRB*RCL_X1S*RCL_CS*ZGAMMA2S
RCL_CONST8S = (3.0_JPRB+RCL_DS+RCL_X4S-RCL_X2S)/(RCL_BS+1.0_JPRB+RCL_X4S-RCL_X2S)

!--------------------------------------------------------
! Define precipitation type
!--------------------------------------------------------
! WMO code table (4.201) for precipitation type
! 0 Reserved
! 1 Rain
! 2 Thunderstorm
! 3 Freezing Rain
! 4 Mixed/Ice
! 5 Snow
! 6 Wet snow
! 7 Melting snow (sleet)
! 8 Ice pellets
! 9 Graupel
! 10 Hail
! 11 Drizzle
! 12 Freezing drizzle

! NPTYPE_SEV2WMO is dimensioned with NPRECTYPES in yoecldp

! Look up table to convert precip type severity to WMO code
! In order of decreasing severity
NPTYPE_SEV2WMO(1)  = 11 ! Drizzle
NPTYPE_SEV2WMO(2)  = 1  ! Rain
NPTYPE_SEV2WMO(3)  = 2  ! Thunderstorm (if class as very heavy rain)
NPTYPE_SEV2WMO(4)  = 7  ! Melting snow (sleet)
NPTYPE_SEV2WMO(5)  = 4  ! Mixed/Ice
NPTYPE_SEV2WMO(6)  = 8  ! Ice pellets
NPTYPE_SEV2WMO(7)  = 9  ! Graupel
NPTYPE_SEV2WMO(8)  = 5  ! Dry snow
NPTYPE_SEV2WMO(9)  = 6  ! Wet snow
NPTYPE_SEV2WMO(10) = 10 ! Hail
NPTYPE_SEV2WMO(11) = 12 ! Freezing drizzle
NPTYPE_SEV2WMO(12) = 3  ! Freezing rain

!----------------------------------------------------------
! Control Cloudsc diagnostics
! Add up these options to get the integer for NCLDDIAG
! 0=no diagnostics
! 1=place microphysics pathways into PEXTRA extra fields
! 2=Perform diagnostics for total water and enthalpy
!----------------------------------------------------------
NCLDDIAG=0

! ---------------------------------------------------------------------
! LCLDBUD logicals store enthalpy and cloud water budgets
! Switches currently hardwired here
! LCLDBUDC        - True = Turn on 3D cloud fraction process budget
! LCLDBUDL        - True = Turn on 3D cloud liquid process budget    
! LCLDBUDI        - True = Turn on 3D cloud ice process budget   
! LCLDBUDT        - True = Turn on 3D cloud process temperature budget   
! LCLDBUD_VERTINT - True = Turn on vertical integrated budget for all terms
! LCLDBUD_TIMEINT - True = Accumulate budget rather than instantaneous. 
!                           Applies to all terms above.
! ---------------------------------------------------------------------
! Default to false. All diagnostics off, can be overridden by namelist.
LCLDBUDC        = .FALSE.
LCLDBUDL        = .FALSE.
LCLDBUDI        = .FALSE.
LCLDBUDT        = .FALSE.
LCLDBUD_VERTINT = .FALSE.
LCLDBUD_TIMEINT = .TRUE.

!-------------------------------
! Aerosol-cloud indirect effects
!-------------------------------
NAERCLD=0 ! 0 = no aerosol interactions
NAECLBC=9
NAECLDU=4
NAECLOM=7
NAECLSS=1
NAECLSU=11
LAERLIQAUTOLSP = .FALSE.
LAERLIQAUTOCP  = .FALSE.
LAERLIQAUTOCPB = .FALSE.
LAERLIQCOLL    = .FALSE.
LAERICESED     = .FALSE.
LAERICEAUTO    = .FALSE.

!-- relationship aerosol mass - CCN numbers (Menon et al., 2002, JAS 59, 695)
RCCNOM = 0.13_JPRB
RCCNSS = 0.05_JPRB
RCCNSU = 0.50_JPRB 

!---------------------
! variance definitions
!---------------------
NSHAPEP=1.0_JPRB+SQRT(2.0_JPRB)
NSHAPEQ=1.0_JPRB+SQRT(2.0_JPRB)
NBETA=100 ! must equal array size of RBETA
DO IX=0,NBETA
  ZX=REAL(IX)/NBETA
  ZX=MAX(MIN(ZX,1.0_JPRB),0.0_JPRB)
  YDECLDP%RBETA(IX)  =BETAI(NSHAPEP,NSHAPEQ,ZX)
  YDECLDP%RBETAP1(IX)=BETAI(NSHAPEP+1,NSHAPEQ,ZX)
ENDDO

! Calculate model level above which cloud scheme is not called
! using a "standard atmosphere" profile
NCLDTOP=2
DO JLEV=NFLEVG,2,-1
  IF (STPRE(JLEV) > RCLDTOPP) NCLDTOP=JLEV
ENDDO

!------------------------------------------------------------
! override the options for the cloud scheme from the namelist
!------------------------------------------------------------
CALL POSNAM(NULNAM,'NAMCLDP')
READ (NULNAM,NAMCLDP)

!--------------------------------------------
! diagnostic switches based on NCLDDIAG
! 1=Extra fields defined 
! 2=Perforce Enthalpy and total water budgets
!--------------------------------------------
LCLDEXTRA  =BTEST(NCLDDIAG,0)
LCLDBUDGET =BTEST(NCLDDIAG,1)

!----------------------------------------------------
! Timescale for supersaturation to cause overcast sky
! Limit to a minimum timescale of one timestep
! Calculated after namelist read as RKOOPTAU might change
!----------------------------------------------------
RSSICEFACTOR = MIN(TSPHY/RKOOPTAU,1.0_JPRB)

!------------------------
! output the final values
!------------------------
WRITE(NULOUT,'("SUCLDP: RCLDIFF =",E12.4)')    RCLDIFF
WRITE(NULOUT,'("SUCLDP: RCLCRIT =",E12.4)')    RCLCRIT
WRITE(NULOUT,'("SUCLDP: RCLCRIT_SEA =",E12.4)')  RCLCRIT_SEA
WRITE(NULOUT,'("SUCLDP: RCLCRIT_LAND =",E12.4)') RCLCRIT_LAND
WRITE(NULOUT,'("SUCLDP: RKCONV =",E12.4)')     RKCONV
WRITE(NULOUT,'("SUCLDP: RLCRITSNOW =",E12.4)') RLCRITSNOW
WRITE(NULOUT,'("SUCLDP: RSNOWLIN1 =",E12.4)')  RSNOWLIN1
WRITE(NULOUT,'("SUCLDP: RSNOWLIN2 =",E12.4)')  RSNOWLIN2
WRITE(NULOUT,'("SUCLDP: RVICE =",E12.4)')      RVICE
WRITE(NULOUT,'("SUCLDP: RVRAIN =",E12.4)')     RVRAIN
WRITE(NULOUT,'("SUCLDP: RVSNOW =",E12.4)')     RVSNOW
WRITE(NULOUT,'("SUCLDP: RICEHI1 =",E12.4)')    RICEHI1
WRITE(NULOUT,'("SUCLDP: RICEHI2 =",E12.4)')    RICEHI2
WRITE(NULOUT,'("SUCLDP: RICEINIT =",E12.4)')   RICEINIT
WRITE(NULOUT,'("SUCLDP: RTHOMO =",E12.4)')     RTHOMO
WRITE(NULOUT,'("SUCLDP: RCOVPMIN =",E12.4)')   RCOVPMIN
WRITE(NULOUT,'("SUCLDP: RPECONS =",E12.4)')    RPECONS
WRITE(NULOUT,'("SUCLDP: RTAUMEL =",E12.4)')    RTAUMEL
WRITE(NULOUT,'("SUCLDP: RKOOPTAU =",E12.4)')   RKOOPTAU
WRITE(NULOUT,'("SUCLDP: RSSICEFACTOR =",E12.4)') RSSICEFACTOR
WRITE(NULOUT,'("SUCLDP: RPRECRHMAX =",E12.4)') RPRECRHMAX
WRITE(NULOUT,'("SUCLDP: RCLDTOPCF =",E12.4)')  RCLDTOPCF    
WRITE(NULOUT,'("SUCLDP: RDEPLIQREFRATE =",E12.4)')  RDEPLIQREFRATE
WRITE(NULOUT,'("SUCLDP: RDEPLIQREFDEPTH =",E12.4)') RDEPLIQREFDEPTH
WRITE(NULOUT,'("SUCLDP: RTAU_CLD_TLAD =",E12.4)') RTAU_CLD_TLAD
WRITE(NULOUT,'("SUCLDP: NSSOPT =",I6)')   NSSOPT
WRITE(NULOUT,'("SUCLDP: NCLDDIAG =",I6)') NCLDDIAG
IF (LCLDEXTRA) WRITE(NULOUT,'("SUCLDP: LCLDEXTRA IS TRUE")') 
IF (LCLDBUDGET) WRITE(NULOUT,'("SUCLDP: LCLDBUDGET IS TRUE")') 
IF (LCLOUD_INHOMOG) WRITE(NULOUT,'("SUCLDP: LCLOUD_INHOMOG IS TRUE")') 
IF (LCLDBUDC) WRITE(NULOUT,'("SUCLDP: LCLDBUDC IS TRUE")')      
IF (LCLDBUDL) WRITE(NULOUT,'("SUCLDP: LCLDBUDL IS TRUE")')      
IF (LCLDBUDI) WRITE(NULOUT,'("SUCLDP: LCLDBUDI IS TRUE")')      
IF (LCLDBUDT) WRITE(NULOUT,'("SUCLDP: LCLDBUDT IS TRUE")')      
IF (LCLDBUD_VERTINT) WRITE(NULOUT,'("SUCLDP: LCLDBUD_VERTINT IS TRUE")')
IF (LCLDBUD_TIMEINT) WRITE(NULOUT,'("SUCLDP: LCLDBUD_TIMEINT IS TRUE")')

LAERLIQAUTOLSP = BTEST(NAERCLD,0)
LAERLIQAUTOCP  = BTEST(NAERCLD,1)
LAERLIQAUTOCPB = BTEST(NAERCLD,2)
LAERLIQCOLL    = BTEST(NAERCLD,3)
LAERICESED     = BTEST(NAERCLD,4)
LAERICEAUTO    = BTEST(NAERCLD,5)

! check: Can't not have both LAERLIQAUTOCP and LAERLIQAUTOCPB
IF (LAERLIQAUTOCPB .AND. LAERLIQAUTOCP) THEN
  WRITE(NULOUT,'("Resetting LAERLIQAUTOCP, choose one convective option only!")')
  LAERLIQAUTOCP=.FALSE.
ENDIF

WRITE(NULOUT,'("SUCLDP: NAERCLD=",I2)')NAERCLD
IF (LAERLIQAUTOLSP) WRITE(NULOUT,'("SUCLDP: LAERLIQAUTOLSP ON")')
IF (LAERLIQAUTOCP)  WRITE(NULOUT,'("SUCLDP: LAERLIQAUTOCP  ON")') 
IF (LAERLIQAUTOCPB) WRITE(NULOUT,'("SUCLDP: LAERLIQAUTOCPB ON")') 
IF (LAERLIQCOLL)    WRITE(NULOUT,'("SUCLDP: LAERLIQCOLL    ON")') 
IF (LAERICESED)     WRITE(NULOUT,'("SUCLDP: LAERICESED     ON")') 
IF (LAERICEAUTO)    WRITE(NULOUT,'("SUCLDP: LAERICEAUTO    ON")') 
WRITE(NULOUT,'("SUCLDP: NCCNOM=",F8.3,2X,"NCCNSS=",F8.3,2X,"NCCNSU=",F8.3)') RCCNOM,RCCNSS,RCCNSU
WRITE(NULOUT,'("SUCLDP: NAECLBC=",I2,1X,"NAECLDU=",I2,1X,"NAECLOM=",I2,1X,"NAECLSS=",I2,1X,&
  & "NAECLSU=",I2)') NAECLBC,NAECLDU,NAECLOM,NAECLSS,NAECLSU

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCLDP',1,ZHOOK_HANDLE)

!-----------------------------------------------------------------------
CONTAINS
FUNCTION BETAI(PP,PQ,PX)

  ! Description:
  !
  ! USES betacf,gammln Returns the incomplete beta function I x (a; b).  
  !
  ! Method:
  !   See Numerical Recipes (Fortran)
  !
  ! Author:
  !   A. Tompkins 2000
  !
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE

  ! input parameters:
REAL(KIND=JPRB)   ,INTENT(IN)    :: PP ! beta shape parameters
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ ! beta shape parameters
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX ! integration limit  

  !  local scalars: 
REAL(KIND=JPRB) :: Z_BT, BETAI 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUCLDP:BETAI',0,ZHOOK_HANDLE)

IF(PX>0.0_JPRB .AND. PX<1.0_JPRB )THEN !FACTORS IN FRONT OF THE CONTINUED FRACTION. 
  Z_BT=EXP(GAMMLN(PP+PQ)-GAMMLN(PP)-GAMMLN(PQ)+PP*LOG(PX)+PQ*LOG(1.0_JPRB-PX)) 
ELSE  
  Z_BT=0.0_JPRB 
ENDIF 
IF(PX<(PP+1.0_JPRB)/(PP+PQ+2.0_JPRB))THEN ! USE CONTINUED FRACTION DIRECTLY.
  BETAI=Z_BT*BETACF(PP,PQ,PX)/PP 
ELSE 
  BETAI=1.0_JPRB-Z_BT*BETACF(PQ,PP,1.0_JPRB-PX)/PQ !
!   use continued fraction after making the symmetry transformation. 
ENDIF 

IF (LHOOK) CALL DR_HOOK('SUCLDP:BETAI',1,ZHOOK_HANDLE)

END FUNCTION BETAI

!-----------------------------------------------------------------------

FUNCTION BETACF(PP,PQ,PX) 
  ! Description:
  !
  !  used by betai: evaluates continued fraction for incomplete 
  !  beta function by modified lentz's method ( x 5.2). 
  !  first step of lentz's method. 
  !
  ! Method:
  !   See Numerical Recipes (Fortran)
  !
  ! Author:
  !   A. Tompkins 2000
  !
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE

  ! input parameters:
REAL(KIND=JPRB)   ,INTENT(IN)    :: PP  ! beta shape parameters
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ  ! beta shape parameters
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX  ! integration limit  

INTEGER(KIND=JPIM) :: I_MAXIT, I_M,I_M2 
REAL(KIND=JPRB) :: ZEPS,Z_FPMIN,Z_AA,Z_C,Z_D,Z_DEL,Z_H,Z_QAB,Z_QAM,Z_QAP, BETACF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUCLDP:BETACF',0,ZHOOK_HANDLE)

I_MAXIT=100
ZEPS = 3.E-7_JPRB
Z_FPMIN=1.E-30_JPRB

Z_QAB=PP+PQ 

!  these q's will be used in factors that occur in the coe cients (6.4.6). 

Z_QAP=PP+1.0_JPRB
Z_QAM=PP-1.0_JPRB 
Z_C=1.0_JPRB 
Z_D=1.0_JPRB-Z_QAB*PX/Z_QAP 
IF(ABS(Z_D)<Z_FPMIN)Z_D=Z_FPMIN 
Z_D=1.0_JPRB/Z_D 
Z_H=Z_D 
I_M=1.0_JPRB
Z_DEL = 2.0_JPRB 
DO WHILE (ABS(Z_DEL-1.0_JPRB)>ZEPS)
  I_M2=2*I_M 
  Z_AA=I_M*(PQ-I_M)*PX/((Z_QAM+I_M2)*(PP+I_M2))
  Z_D=1.0_JPRB+Z_AA*Z_D  ! one step (the even one) of the recurrence. 
  IF (ABS(Z_D)<Z_FPMIN)Z_D=Z_FPMIN 
  Z_C=1.0_JPRB+Z_AA/Z_C 
  IF (ABS(Z_C)<Z_FPMIN)Z_C=Z_FPMIN 
  Z_D=1.0_JPRB/Z_D 
  Z_H=Z_H*Z_D*Z_C 
  Z_AA=-(PP+I_M)*(Z_QAB+I_M)*PX/((PP+I_M2)*(Z_QAP+I_M2)) 
  Z_D=1.0_JPRB+Z_AA*Z_D ! next step of the recurrence (the odd one). 
  IF (ABS(Z_D)<Z_FPMIN)Z_D=Z_FPMIN 
  Z_C=1.0_JPRB+Z_AA/Z_C 
  IF(ABS(Z_C)<Z_FPMIN)Z_C=Z_FPMIN 
  Z_D=1.0_JPRB/Z_D 
  Z_DEL=Z_D*Z_C 
  Z_H=Z_H*Z_DEL 
  I_M=I_M+1.0_JPRB !*AMT*
  IF(I_M>I_MAXIT) THEN 
    Z_DEL=1.0_JPRB
  ENDIF  
ENDDO 
BETACF=Z_H 

IF (LHOOK) CALL DR_HOOK('SUCLDP:BETACF',1,ZHOOK_HANDLE)
END FUNCTION BETACF

!-----------------------------------------------------------------------

FUNCTION GAMMLN(P_XX) 

  ! Description:
  !
  ! Gamma function calculation
  ! returns the value ln[g(xx)] for xx > _ZERO_
  !
  ! Method:
  !   See Numerical Recipes
  !
  ! Author:
  !   A. Tompkins 2000
  !
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: P_XX

REAL(KIND=JPRB) :: GAMMLN
INTEGER(KIND=JPIM) :: J 
REAL(KIND=JPRB) :: Z_SER,Z_STP,Z_TMP,Z_X,Z_Y,Z_COF(6) 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUCLDP:GAMMLN',0,ZHOOK_HANDLE)

Z_COF(:)= (/&
 & 76.18009172947146_JPRB, &
 & -86.50532032941677_JPRB, &
 & 24.01409824083091_JPRB, &
 & -1.231739572450155_JPRB, &
 & 0.1208650973866179E-2_JPRB, &
 & -.5395239384953E-5_JPRB/)  

Z_STP=2.5066282746310005_JPRB

Z_X=P_XX 
Z_Y=Z_X 
Z_TMP=Z_X+5.5_JPRB
Z_TMP=(Z_X+0.5_JPRB)*LOG(Z_TMP)-Z_TMP 
Z_SER=1.000000000190015_JPRB
DO J=1,6
  Z_Y=Z_Y+1.0_JPRB 
  Z_SER=Z_SER+Z_COF(J)/Z_Y 
ENDDO 
GAMMLN=Z_TMP+LOG(Z_STP*Z_SER/Z_X) 

IF (LHOOK) CALL DR_HOOK('SUCLDP:GAMMLN',1,ZHOOK_HANDLE)
END FUNCTION GAMMLN
END SUBROUTINE SUCLDP
