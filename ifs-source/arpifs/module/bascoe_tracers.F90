! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE BASCOE_TRACERS

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

  INTEGER(KIND=JPIM), PARAMETER  :: ntrace = 64  ! All tracers for chemistry
  !
  ! components numbers
  !
  INTEGER(KIND=JPIM), PARAMETER :: iO3      =1
  INTEGER(KIND=JPIM), PARAMETER :: iH2O2    =2 
  INTEGER(KIND=JPIM), PARAMETER :: iCH4     =3
  INTEGER(KIND=JPIM), PARAMETER :: iCO      =4
  INTEGER(KIND=JPIM), PARAMETER :: iHNO3    =5
  INTEGER(KIND=JPIM), PARAMETER :: iCH3OOH  =6 
  INTEGER(KIND=JPIM), PARAMETER :: iCH2O    =7
  INTEGER(KIND=JPIM), PARAMETER :: iNO      =8
  INTEGER(KIND=JPIM), PARAMETER :: iHO2     =9
  INTEGER(KIND=JPIM), PARAMETER :: iCH3    =10
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O   =11
  INTEGER(KIND=JPIM), PARAMETER :: iHCO    =12
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O2  =13 
  INTEGER(KIND=JPIM), PARAMETER :: iOH     =14
  INTEGER(KIND=JPIM), PARAMETER :: iNO2    =15
  INTEGER(KIND=JPIM), PARAMETER :: iN2O5   =16
  INTEGER(KIND=JPIM), PARAMETER :: iHO2NO2 =17
  INTEGER(KIND=JPIM), PARAMETER :: iNO3    =18
  INTEGER(KIND=JPIM), PARAMETER :: iN2O    =19
  INTEGER(KIND=JPIM), PARAMETER :: iH2O    =20
  INTEGER(KIND=JPIM), PARAMETER :: iOCLO   =21
  INTEGER(KIND=JPIM), PARAMETER :: iHCL    =22
  INTEGER(KIND=JPIM), PARAMETER :: iCLONO2 =23
  INTEGER(KIND=JPIM), PARAMETER :: iHOCL   =24
  INTEGER(KIND=JPIM), PARAMETER :: iCL2    =25
  INTEGER(KIND=JPIM), PARAMETER :: iHBR    =26
  INTEGER(KIND=JPIM), PARAMETER :: iBRONO2 =27 
  INTEGER(KIND=JPIM), PARAMETER :: iCL2O2  =28
  INTEGER(KIND=JPIM), PARAMETER :: iHOBR   =29
  INTEGER(KIND=JPIM), PARAMETER :: iBRCL   =30
  INTEGER(KIND=JPIM), PARAMETER :: iCFC11  =31
  INTEGER(KIND=JPIM), PARAMETER :: iCFC12  =32
  INTEGER(KIND=JPIM), PARAMETER :: iCFC113 =33
  INTEGER(KIND=JPIM), PARAMETER :: iCFC114 =34
  INTEGER(KIND=JPIM), PARAMETER :: iCFC115 =35
  INTEGER(KIND=JPIM), PARAMETER :: iCCL4   =36
  INTEGER(KIND=JPIM), PARAMETER :: iCLNO2  =37
  INTEGER(KIND=JPIM), PARAMETER :: iCH3CCL3=38
  INTEGER(KIND=JPIM), PARAMETER :: iCH3CL  =39
  INTEGER(KIND=JPIM), PARAMETER :: iHCFC22 =40
  INTEGER(KIND=JPIM), PARAMETER :: iCH3BR  =41
  INTEGER(KIND=JPIM), PARAMETER :: iHF     =42
  INTEGER(KIND=JPIM), PARAMETER :: iHA1301 =43
  INTEGER(KIND=JPIM), PARAMETER :: iHA1211 =44
  INTEGER(KIND=JPIM), PARAMETER :: iCHBR3  =45
  INTEGER(KIND=JPIM), PARAMETER :: iCLOO   =46
  INTEGER(KIND=JPIM), PARAMETER :: iO      =47
  INTEGER(KIND=JPIM), PARAMETER :: iO1D    =48
  INTEGER(KIND=JPIM), PARAMETER :: iN      =49
  INTEGER(KIND=JPIM), PARAMETER :: iCLO    =50
  INTEGER(KIND=JPIM), PARAMETER :: iCL     =51
  INTEGER(KIND=JPIM), PARAMETER :: iBR     =52
  INTEGER(KIND=JPIM), PARAMETER :: iBRO    =53
  INTEGER(KIND=JPIM), PARAMETER :: iH      =54
  INTEGER(KIND=JPIM), PARAMETER :: iH2     =55
  INTEGER(KIND=JPIM), PARAMETER :: iCO2    =56
  INTEGER(KIND=JPIM), PARAMETER :: iBR2    =57
  INTEGER(KIND=JPIM), PARAMETER :: ICH2BR2 =58
  INTEGER(KIND=JPIM), PARAMETER :: iStratAer=59
  INTEGER(KIND=JPIM), PARAMETER :: IOCS    =60
  INTEGER(KIND=JPIM), PARAMETER :: ISO2    =61
  INTEGER(KIND=JPIM), PARAMETER :: ISO3    =62
  INTEGER(KIND=JPIM), PARAMETER :: ISO4    =63
  INTEGER(KIND=JPIM), PARAMETER :: IH2SO4  =64


  ! Boundary conditions at surface for stratospheric species for IFS(BASCOE)-only.
  ! Which currently don't feature emissions (!) 
  ! The list used for IFS(CB05BASCOE) experiments is in 'bascoetm5_module'
  ! Note that we also introduce BC for O3 and CH4, to prevent a blow-up in the troposphere
  INTEGER(KIND=JPIM), PARAMETER  :: NBC = 21
  INTEGER(KIND=JPIM),DIMENSION(NBC), PARAMETER ::BASCOE_BC=(/ &
  & IN2O  , ICFC11  , ICFC12 , ICFC113, ICFC114, &
  & ICFC115, ICH2BR2, ICCL4 , ICH3CCL3, IHCFC22, &
  & IHA1301, IHA1211, ICH3Br, ICHBR3  , ICH3CL , &
  & ICO2   , IO3    , ICH4  , IBRO    , IHBR   , &
  & IOCS /)
  CHARACTER (LEN = 12), DIMENSION(NBC) :: BASCOE_BCNAME = &
     &      (/"n2o         " , &
     &        "cfc11       " , &
     &        "cfc12       " , &
     &        "cfc113      " , &
     &        "cfc114      " , &
     &        "cfc115      " , &
     &        "ch2br2      " , &
     &        "ccl4        " , &
     &        "ch3ccl3     " , &
     &        "hcfc22      " , &
     &        "ha1301      " , &
     &        "ha1211      " , &
     &        "ch3br       " , &
     &        "chbr3       " , &
     &        "ch3cl       " , &
     &        "co2         " , &
     &        "o3          " , &
     &        "ch4         " , &
     &        "bro         " , &
     &        "hbr         " , &
     &        "ocs         " /)

  ! Just a simple parameter, taken from BASCOE field for 1 April 2008
  ! Units in [mole/mole]
  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL=(/ &
  &  3.27E-7 , 2.33E-10, 5.20E-10, 7.27E-11, 1.63E-11, &
  &  8.43E-12, 1.20E-12, 83.2E-12, 3.70E-12, 2.29E-10, &
  &  3.30E-12, 3.75E-12, 6.65E-12, 1.17E-12, 5.34E-10, &
  &  400.E-6 , 3.00E-8 , 1.74E-6 , 0.10E-12, 0.10E-12, &
  &  266.E-12 /)
  ! *WAERNING* o3 here is arbitrary !
  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL_1750=(/ &
  &  2.74E-7 , 0.00E-10, 0.00E-10, 0.00E-11, 0.00E-12, &
  &  0.00E-12, 1.20E-12, 2.94E-14, 0.00E-11, 0.00E-10, &
  &  0.00E-12, 0.00E-12, 5.30E-12, 1.20E-12, 4.57E-10, &
  &  277.E-6 , 3.00E-8 , 0.731E-6, 0.00E-12, 0.00E-12, &
  &  266.E-12 /)

  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL_1850=(/ &
  &  2.73E-7 , 0.00E-10, 0.00E-10, 0.00E-11, 0.00E-12, &
  &  0.00E-10, 1.20E-12, 2.94E-14, 0.00E-11, 0.00E-10, &
  &  0.00E-12, 0.00E-12, 5.30E-12, 1.20E-12, 4.57E-10, &
  &  284.E-6 , 3.00E-8 , 0.808E-6, 0.00E-12, 0.00E-12, &
  &  266.E-12 /)


END MODULE BASCOE_TRACERS
