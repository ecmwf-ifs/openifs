! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_AERP(YDEAERATM,YDML_PHY_AER,YDCOMPO)

!**** *SU_AERP*   - INITIALIZE MODULES YOEAERSRC, YOEAERSNK

!     PURPOSE.
!     --------
!           INITIALIZE YOEAERSRC AND YOEAERSNK, THE MODULES THAT CONTAINS 
!           COEFFICIENTS NEEDED TO RUN THE PROGNOSTIC AEROSOLS

!**   INTERFACE.
!     ----------
!        *CALL* *SU_AERP

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        YOEAERSRC, YOEAERSNK, YOEAERATM

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE *ECMWF*
!        from O.BOUCHER (LOA, 1998-03) 

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2004-05-10
!        JJMorcrette 20070222 9bins for SS and DU for testing
!        SRémy       20160308 Update of conversion rates for SO2 to SO4
!        SRémy       20160621 dry matter SS, update of RHO, RMMD
!        SRémy       20160916 Add nitrates
!        SRémy       20170420 Wet deposition fraction for SS increased (0.2 -> 0.9)
!        SRémy       20171113 Add SOA
!        SRémy       20230403 Use Gong03 for SS aerosol emissions
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_AEROSOL_MOD , ONLY : MODEL_PHYSICS_AEROSOL_TYPE
USE PARKIND1                  , ONLY : JPRB , JPIM
USE YOMHOOK                   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOEAEROP                  ,ONLY :  CORRALF_DD, CORRASY_DD, CORROMG_DD



USE YOEAERATM , ONLY : TEAERATM
USE YOMCOMPO  , ONLY : TCOMPO

IMPLICIT NONE


!*       0.5   LOCAL VARIABLES
!              ---------------

TYPE(TEAERATM)                  ,INTENT(INOUT) :: YDEAERATM
TYPE(MODEL_PHYSICS_AEROSOL_TYPE),INTENT(INOUT) :: YDML_PHY_AER
TYPE(TCOMPO)                    ,INTENT(INOUT) :: YDCOMPO
INTEGER(KIND=JPIM) :: IDU, ISS

!-- sea-salt (3 or 9 bins)
REAL(KIND=JPRB) :: ZSSFLX3(3)
REAL(KIND=JPRB) :: ZSSFLX9(9)
REAL(KIND=JPRB) :: ZMMD_SS3(3)
REAL(KIND=JPRB) :: ZMMD_SS9(9)

!-- desert dust (3 or 9 bins)
REAL(KIND=JPRB) :: ZMMD_DD3(3),ZRHO_DD3(3)
REAL(KIND=JPRB) :: ZMMD_DD9(9),ZRHO_DD9(9)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_AERP',0,ZHOOK_HANDLE)
ASSOCIATE(RGRATE=>YDEAERATM%RGRATE, &
 & LAERNITRATE=>YDCOMPO%LAERNITRATE, &
 & LAERSOA=>YDCOMPO%LAERSOA, &
 & LSEASALT_RH80=>YDEAERATM%LSEASALT_RH80, &
 & LAERDUST_NEWBIN=>YDEAERATM%LAERDUST_NEWBIN, &
 & RFRAER=>YDML_PHY_AER%YREAERSNK%RFRAER, &
 & RMMD_STRATSULF_BACK=>YDML_PHY_AER%YREAERSNK%RMMD_STRATSULF_BACK, &
 & RMMD_STRATSULF_VOLC=>YDML_PHY_AER%YREAERSNK%RMMD_STRATSULF_VOLC, &
 & RSIGMA_STRATSULF_VOLC=>YDML_PHY_AER%YREAERSNK%RSIGMA_STRATSULF_VOLC, &
 & RRHO_STRATSULF=>YDML_PHY_AER%YREAERSNK%RRHO_STRATSULF, &
 & RHO_ICE=>YDML_PHY_AER%YREAERSNK%RHO_ICE, RHO_WAT=>YDML_PHY_AER%YREAERSNK%RHO_WAT, &
 & RMMD_DD=>YDML_PHY_AER%YREAERSNK%RMMD_DD, &
 & RMMD_NI=>YDML_PHY_AER%YREAERSNK%RMMD_NI, &
 & RMMD_AM=>YDML_PHY_AER%YREAERSNK%RMMD_AM, &
 & RMMD_SS=>YDML_PHY_AER%YREAERSNK%RMMD_SS, &
 & RMMD_SOA=>YDML_PHY_AER%YREAERSNK%RMMD_SOA, &
 & RMMD_SO4=>YDML_PHY_AER%YREAERSNK%RMMD_SO4, &
 & RMMD_OM=>YDML_PHY_AER%YREAERSNK%RMMD_OM, &
 & RMMD_BC=>YDML_PHY_AER%YREAERSNK%RMMD_BC, &
 & RMMD_ASH=>YDML_PHY_AER%YREAERSNK%RMMD_ASH, &
 & RRHO_DD=>YDML_PHY_AER%YREAERSNK%RRHO_DD, &
 & RRHO_NI=>YDML_PHY_AER%YREAERSNK%RRHO_NI, &
 & RRHO_AM=>YDML_PHY_AER%YREAERSNK%RRHO_AM, &
 & RRHO_SS=>YDML_PHY_AER%YREAERSNK%RRHO_SS, &
 & RRHO_SOA=>YDML_PHY_AER%YREAERSNK%RRHO_SOA, &
 & RRHO_SO4=>YDML_PHY_AER%YREAERSNK%RRHO_SO4, &
 & RRHO_OM=>YDML_PHY_AER%YREAERSNK%RRHO_OM, &
 & RRHO_BC=>YDML_PHY_AER%YREAERSNK%RRHO_BC, &
 & RRHO_ASH=>YDML_PHY_AER%YREAERSNK%RRHO_ASH, &
 & RRHTAB=>YDML_PHY_AER%YREAERSNK%RRHTAB, RSSGROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%RSSGROWTH_RHTAB, &
 & ROMGROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%ROMGROWTH_RHTAB, &
 & RSOAGROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%RSOAGROWTH_RHTAB, &
 & RNIGROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%RNIGROWTH_RHTAB, &
 & RAMGROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%RAMGROWTH_RHTAB, &
 & RSO4GROWTH_RHTAB=>YDML_PHY_AER%YREAERSNK%RSO4GROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YDML_PHY_AER%YREAERSNK%RSSDENS_RHTAB, &
 & RHAMAKER=>YDML_PHY_AER%YREAERSNK%RHAMAKER, &
 & RSO2CV1=>YDML_PHY_AER%YREAERSNK%RSO2CV1, &
 & RSO2CV2=>YDML_PHY_AER%YREAERSNK%RSO2CV2, RSUCV1=>YDML_PHY_AER%YREAERSNK%RSUCV1, &
 & RSUCV2=>YDML_PHY_AER%YREAERSNK%RSUCV2, &
 & RVSO2CV1=>YDML_PHY_AER%YREAERSNK%RVSO2CV1, RVSO2CV2=>YDML_PHY_AER%YREAERSNK%RVSO2CV2, &
 & R_R=>YDML_PHY_AER%YREAERSNK%R_R, &
 & R_S=>YDML_PHY_AER%YREAERSNK%R_S, &
 & NTYPAER=>YDEAERATM%NTYPAER, &
 & RSSFLX=>YDML_PHY_AER%YREAERSRC%RSSFLX)

!-- For the ECMWF model, the following tables when dimensioned to 12 
!   can refer to 12 values of RH
!      (RHTAB, RHHO_SS)
!   or to 12 types/bins of aerosols with the following mapping: NTYPAER
!     1       1- 3  sea-salt  0.03 - 0.5 -  5  - 20 microns
!     2       4- 6  dust      0.03 - 0.5 - 0.9 - 20 microns
!     3       7- 8  POM    hydrophilic, hydrophobic
!     4       9-10  BC     hydrophilic, hydrophobic
!     5      11-12  SO4/SO2 including sulfate prognostic stratospheric aerosols
!     (SO4 is 11)
!     6      13-14  Nitrate fine and coarse
!     7      15     Ammonium
!     8      16     fly ash
!     9      17-18   volcanic SO2/SO4
!      (RVDPOCE, RVDSIC, RVDPLND, RVDPLIC)
!      (RVSEDOCE)

!-- parameters are:
! RFRxx       efficiency for in-cloud scavenging
! RALPHAR, S  efficiency for below-cloud scavenging
! RVDPyy      speed relevant to dry deposition
! RVSEDzz     speed relevant to sedimentation



!*      0.    INITIALISATION
!             --------------

RSSFLX(:)  = 0._JPRB

RRHO_DD(:) = 0._JPRB
RRHO_SS(:) = 0._JPRB
RRHO_NI(:) = 0._JPRB

RRHO_OM = 0._JPRB
RRHO_AM = 0._JPRB
RRHO_BC = 0._JPRB
RRHO_SOA = 0._JPRB
RRHO_SO4 = 0._JPRB
RRHO_ASH = 0._JPRB

RMMD_NI(:) = (/0._JPRB,0._JPRB/)
RMMD_OM = 0._JPRB
RMMD_AM = 0._JPRB
RMMD_BC = 0._JPRB
RMMD_SOA = 0._JPRB
RMMD_SO4 = 0._JPRB
RMMD_ASH = 0._JPRB
RMMD_STRATSULF_BACK=0.8_JPRB
RMMD_STRATSULF_VOLC=0.8_JPRB
RSIGMA_STRATSULF_VOLC=1.6_JPRB
RRHO_STRATSULF=1730._JPRB

!*      1.    PARAMETERS RELATED TO AEROSOL TYPES
!             -----------------------------------

! NMAXTAER is now PARAMETER in YOEAERATM


!*      2.    PARAMETERS RELATED TO SOURCES/SINKS
!             -----------------------------------

!- Rules of thummb for potential adjustment:
! *  FRxx efficiency parameter for in-cloud scavenging 
!         the greater, the more efficient (via an negative exponential)
! *  ALPHAR, ALPHAS efficiency for below-cloud scavenging
!         the greater, the more efficient (directly proportional)
! *  VSED sedimentation coefficient
!         the greater, the larger the sedimentation flux
! *  VDEP dry deposition 
!         the greater, the larger the loss in the lowermost layer

R_R = 0.001_JPRB
R_S = 0.001_JPRB

RFRAER = 0.5_JPRB

! Concrete-concrete-air, Kim et al 2010, "Source term models for fine particle resuspension from indoor surface"
RHAMAKER=5.E-20_JPRB

!*      2.1   SEA SALT
!             -------- 
!-- parameters related to SEA SALT: 12 relates to 12 values of relative humidity

RRHTAB = (/ 0._JPRB, 10._JPRB, 20._JPRB, 30._JPRB, 40._JPRB, 50._JPRB &
       & , 60._JPRB, 70._JPRB, 80._JPRB, 85._JPRB, 90._JPRB, 95._JPRB /)

! SS hygroscopic growth (Tang et al 1997)

RSSGROWTH_RHTAB = (/ 1._JPRB, 1._JPRB, 1._JPRB, 1._JPRB, 1.442_JPRB, 1.555_JPRB,  &
                   & 1.666_JPRB, 1.799_JPRB, 1.988_JPRB, 2.131_JPRB, 2.361_JPRB, 2.876_JPRB/)
RSSDENS_RHTAB = (/2160._JPRB ,2160._JPRB ,2160._JPRB,  2160._JPRB,  1451.6_JPRB, 1367.9_JPRB, &
          &       1302.9_JPRB,1243.2_JPRB,1182.7_JPRB, 1149.5_JPRB,1111.6_JPRB, 1063.1_JPRB/)

!* Growth factors corresponding to RH table as given in RRHTAB, according to
!J.J. Morcrette
RSO4GROWTH_RHTAB = &
         & (/1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.169_JPRB,1.220_JPRB,&
         & 1.282_JPRB,1.363_JPRB,1.485_JPRB,1.581_JPRB,1.732_JPRB,2.085_JPRB/)
ROMGROWTH_RHTAB =  &
         & (/1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.05_JPRB,&
         & 1.100_JPRB,1.150_JPRB,1.200_JPRB,1.250_JPRB,1.300_JPRB,1.400_JPRB/)
RNIGROWTH_RHTAB =  &
        & (/1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.100_JPRB,1.200_JPRB,&
        & 1.250_JPRB,1.300_JPRB,1.350_JPRB,1.500_JPRB,1.700_JPRB,2.100_JPRB/)
RAMGROWTH_RHTAB =  &
        & (/1.00_JPRB,1.00_JPRB,1.00_JPRB,1.00_JPRB,1.169_JPRB,1.220_JPRB,&
        & 1.282_JPRB,1.363_JPRB,1.485_JPRB,1.581_JPRB,1.732_JPRB,2.085_JPRB/)
RSOAGROWTH_RHTAB =  &
         & (/1.00_JPRB,1.05_JPRB,1.05_JPRB,1.05_JPRB,1.05_JPRB,1.05_JPRB,&
         & 1.05_JPRB,1.05_JPRB,1.05_JPRB,1.1_JPRB,1.15_JPRB,1.2_JPRB/)

RRHO_AM=1760._JPRB  ! kg/m^3
RRHO_NI=(/1730._JPRB,1400._JPRB/)  ! kg/m^3
RRHO_OM=1800._JPRB  ! kg/m^3
RRHO_BC=1000._JPRB  ! kg/m^3
RRHO_SO4=1760._JPRB  ! kg/m^3
RRHO_SOA=1800._JPRB  ! kg/m^3
RRHO_ASH=2700._JPRB  ! kg/m^3

RHO_WAT = 1000._JPRB
RHO_ICE = 917._JPRB

RMMD_NI(:) = (/0.352_JPRB,2.58_JPRB/)
RMMD_OM = 0.312_JPRB
RMMD_AM = 0.352_JPRB
RMMD_BC = 0.288_JPRB
RMMD_SOA = 0.312_JPRB
RMMD_SO4 = 0.352_JPRB
RMMD_ASH = 1.41_JPRB

!-- PARAMETERS DEPENDENT ON NUMBER OF BINS
!   -------------------------------------- 

!N.B. Fluxes of sea salt for each size bin are given in mg m-2 s-1 at wind 
!     speed of 1 m s-1 at 10m height (at 80% RH) in OB's seasalt.F
!     RSSFLX also in mg m-2 s-1       
! SR 06/2019: updated value
!ZSSFLX3  = (/ 1.63663536E-08_JPRB, 0.9768556E-06_JPRB, 2.02305813E-06_JPRB /)
!SR 04/2023: Gong03
ZSSFLX3  = (/ 1.21837E-08_JPRB, 0.32256E-06_JPRB, 1.18944E-06_JPRB /)

ZMMD_SS3 = (/ 0.68_JPRB, 2.58_JPRB, 12.92_JPRB /)



!N.B. Fluxes of sea salt for each size bin are given in mg m-2 s-1 at wind 
!     speed of 1 m s-1 at 10m height (at 80% RH) in OB's seasalt.F
!     RSSFLX also in mg m-2 s-1       
ZSSFLX9 = (/ 1.36417419E-10_JPRB, 3.46521922E-10_JPRB, 7.53073937E-10_JPRB, 3.62362051E-09_JPRB, &
  &          1.62356013E-08_JPRB, 8.62133902E-08_JPRB, 3.12909776E-07_JPRB, 2.46328341E-07_JPRB, &
  &          2.58577359E-07_JPRB /)
!-- empirically adjusted
ZMMD_SS9 = (/  0.042_JPRB, 0.081_JPRB, 0.170_JPRB, 0.353_JPRB, 0.792_JPRB, &
  &            1.711_JPRB, 3.543_JPRB, 7.951_JPRB,17.122_JPRB /)


!-- empirically adjusted
ZMMD_SS9 = (/  0.042_JPRB, 0.081_JPRB, 0.170_JPRB, 0.353_JPRB, 0.792_JPRB, &
  &            1.711_JPRB, 3.543_JPRB, 7.951_JPRB,17.122_JPRB /)

ISS=NTYPAER(1)

IF (ISS == 9) THEN
  RSSFLX(1:9)  = ZSSFLX9(1:9)
  RMMD_SS(1:9) = ZMMD_SS9(1:9)
ELSE
  RSSFLX(1:3)  = ZSSFLX3(1:3)
  RMMD_SS(1:3) = ZMMD_SS3(1:3)
ENDIF

IF (LSEASALT_RH80) THEN
  RRHO_SS(1:ISS) = RSSDENS_RHTAB(9)
ELSE
  RRHO_SS(1:ISS) = RSSDENS_RHTAB(1)
  RMMD_SS(1:ISS) = RMMD_SS(1:ISS) / RSSGROWTH_RHTAB(9)
ENDIF


!*      2.2   DESERT DUST
!             ----------- 

!-- ECMWF 3 bins of desert dust

ZMMD_DD3(1:3) = (/ 0.51_JPRB, 1.41_JPRB, 5.06_JPRB /)

!- parameters related to DESERT DUST  (ECMWF 3 bins)
!  bins are 0.03 - 0.55 -  0.9  - 20 microns

ZRHO_DD3(1:3) = (/ 2600._JPRB, 2600._JPRB, 2600._JPRB /)

!-- ECMWF 9 bins of desert dust
!  bins are 0.03 - 0.05 - 0.10 - 0.20 - 0.50 - 1.00 - 2.00 - 5.00 - 10.00 - 20.00 microns
ZMMD_DD9 = (/  0.042_JPRB, 0.081_JPRB, 0.170_JPRB, 0.353_JPRB, 0.792_JPRB, &
  &            1.711_JPRB, 3.543_JPRB, 7.951_JPRB,17.122_JPRB /)
ZRHO_DD9 = (/  2600._JPRB, 2600._JPRB, 2600._JPRB, 2600._JPRB, 2600._JPRB, &
  &            2600._JPRB, 2600._JPRB, 2600._JPRB, 2600._JPRB /)


IDU=NTYPAER(2)

IF (IDU == 9) THEN
  RRHO_DD(1:9) = ZRHO_DD9(1:9)
  RMMD_DD(1:9) = ZMMD_DD9(1:9)
ELSE
  RRHO_DD(1:3) = ZRHO_DD3(1:3)
  RMMD_DD(1:3) = ZMMD_DD3(1:3)
ENDIF

! Asphericity correct

CORRALF_DD(:,1) =  (/1.034_JPRB, 1.044_JPRB, 1.061_JPRB, 1.072_JPRB, 1.081_JPRB, &
                  & 1.075_JPRB, 1.063_JPRB, 1.049_JPRB,1.042_JPRB,1.005_JPRB,0.999_JPRB, &
                  & 0.979_JPRB, 0.976_JPRB, 0.975_JPRB, 0.973_JPRB, 0.972_JPRB,0.968_JPRB, 0.953_JPRB, &
                  & 0.980_JPRB, 1._JPRB/)

CORROMG_DD(:,1) =  (/1.006_JPRB, 1.007_JPRB, 1.009_JPRB, 1.007_JPRB, 1.006_JPRB, &
                  & 1.004_JPRB, 1.002_JPRB, 1.002_JPRB,1.0_JPRB,1.0_JPRB,0.999_JPRB, &
                  & 0.999_JPRB, 0.999_JPRB, 0.999_JPRB, 0.999_JPRB, 0.999_JPRB,0.999_JPRB, 0.999_JPRB, &
                  & 0.999_JPRB, 1._JPRB/)

CORRASY_DD(:,1) =  (/1.019_JPRB, 1.027_JPRB, 1.041_JPRB, 1.049_JPRB, 1.058_JPRB, &
                  & 1.056_JPRB, 1.049_JPRB, 1.043_JPRB,1.037_JPRB,1.021_JPRB,1.019_JPRB, &
                  & 1.009_JPRB, 1.007_JPRB, 1.006_JPRB, 1.006_JPRB, 1.005_JPRB,1.007_JPRB, 0.989_JPRB, &
                  & 0.971_JPRB, 1._JPRB/)

CORRALF_DD(:,2) =  (/1.052_JPRB, 1.051_JPRB, 1.052_JPRB, 1.045_JPRB, 1.018_JPRB, &
                  & 1.004_JPRB, 1.005_JPRB, 1.023_JPRB,1.043_JPRB,1.1699_JPRB,1.184_JPRB, &
                  & 1.152_JPRB, 1.110_JPRB, 1.108_JPRB, 1.037_JPRB, 1.020_JPRB,0.990_JPRB, 0.979_JPRB, &
                  & 0.972_JPRB, 1._JPRB/)

CORROMG_DD(:,2) =  (/1.001_JPRB, 1.0_JPRB, 1.001_JPRB, 1.0_JPRB, 0.997_JPRB, &
                  & 0.996_JPRB, 0.998_JPRB, 1.001_JPRB,1.004_JPRB,1.011_JPRB,1.011_JPRB, &
                  & 1.005_JPRB, 1.004_JPRB, 1.004_JPRB, 1.001_JPRB, 1.001_JPRB,0.999_JPRB, 0.999_JPRB, &
                  & 0.999_JPRB, 1._JPRB/)

CORRASY_DD(:,2) =  (/0.999_JPRB, 0.995_JPRB, 0.990_JPRB, 0.985_JPRB, 0.970_JPRB, &
                  & 0.962_JPRB, 0.965_JPRB, 0.877_JPRB,0.994_JPRB,1.108_JPRB,1.121_JPRB, &
                  & 1.092_JPRB, 1.074_JPRB, 1.072_JPRB, 1.029_JPRB, 1.025_JPRB,1.014_JPRB, 1.011_JPRB, &
                  & 1.009_JPRB, 1._JPRB/)

CORRALF_DD(:,3) =  (/1.055_JPRB, 1.054_JPRB, 1.052_JPRB, 1.051_JPRB, 1.055_JPRB, &
                  & 1.059_JPRB, 1.060_JPRB, 1.059_JPRB,1.058_JPRB,1.046_JPRB,1.044_JPRB, &
                  & 1.046_JPRB, 1.051_JPRB, 1.052_JPRB, 1.065_JPRB, 1.069_JPRB,1.075_JPRB, 1.067_JPRB, &
                  & 1.052_JPRB, 1._JPRB/)

CORROMG_DD(:,3) =  (/1.009_JPRB, 1.008_JPRB, 1.006_JPRB, 1.005_JPRB, 1.003_JPRB, &
                 & 1.003_JPRB, 1.003_JPRB, 1.0_JPRB,0.999_JPRB,0.997_JPRB,0.996_JPRB, &
                 & 0.998_JPRB, 0.999_JPRB, 0.999_JPRB, 1.001_JPRB, 1.002_JPRB,1.003_JPRB, 1.003_JPRB, &
                 & 1.004_JPRB, 1._JPRB/)

CORRASY_DD(:,3) =  (/1.001_JPRB, 1.002_JPRB, 1.002_JPRB, 1.003_JPRB, 1.003_JPRB, &
                  & 1.002_JPRB, 1.002_JPRB, 0.999_JPRB,0.996_JPRB,0.988_JPRB,0.987_JPRB, &
                  & 0.987_JPRB, 0.988_JPRB, 0.989_JPRB, 0.999_JPRB, 1.002_JPRB,1.011_JPRB, 1.021_JPRB, &
                  & 1.024_JPRB, 1._JPRB/)



IF (NTYPAER(1) /= 9 .AND. NTYPAER(2) /= 9) THEN

!*      2.3   OTHER AEROSOLS (to be improved later!)
!             --------------
!- parameters related to other aerosol types

!- sulfate (SO4 / SO2)
!-- coefficients describing the transfer from SO2 to SO4 as function of latitude
! -- GEMS original from LMD model, as in Huneeus et al., 2009, $3.1, eq.3, p.215
!  RSO2CV1 = 8.0_JPRB
!  RSO2CV2 = 5.0_JPRB
! -- adjusted     
!  RSO2CV1 = 7.5_JPRB
!  RSO2CV2 = 7.0_JPRB    
! -- from Météo-France
  RSO2CV1 = 4.0_JPRB
  RSO2CV2 = 3.5_JPRB    

!- volcanic sulfate (SO4 / SO2)  NB: EVERYTHING IS AS OF 20130729 SIMILAR TO TROP.SO2/SO4
!-- coefficients describing the transfer from SO2 to SO4 as function of latitude
! -- GEMS original from LMD model, as in Huneeus et al., 2009, $3.1, eq.3, p.215
! NB: the coefficients for the transfer from volcanic SO2 to SO4 are THE SAME as above
  RVSO2CV1 = 8.0_JPRB
  RVSO2CV2 = 5.0_JPRB

! -- adjusted     
!  RVSO2CV1 = 7.5_JPRB
!  RVSO2CV2 = 7.0_JPRB    



ENDIF


!*      3.    PARAMETERS RELATED TO TRANSPORT WITHIN THE FREE ATMOSPHERE
!             ----------------------------------------------------------

!-- related to conversion between -phobic and -philic
RGRATE = 7.1E-05_JPRB

!     ----------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SU_AERP',1,ZHOOK_HANDLE)
END SUBROUTINE SU_AERP

