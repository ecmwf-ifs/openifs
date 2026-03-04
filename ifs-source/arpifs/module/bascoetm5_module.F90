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

MODULE BASCOETM5_MODULE

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

  INTEGER(KIND=JPIM), PARAMETER  :: ntrace = 120  ! All tracers for chemistry

  !
  ! components numbers
  !
  INTEGER(KIND=JPIM), PARAMETER :: iO3        = 1
  INTEGER(KIND=JPIM), PARAMETER :: iH2O2      = 2 
  INTEGER(KIND=JPIM), PARAMETER :: iCH4       = 3
  INTEGER(KIND=JPIM), PARAMETER :: iCO        = 4
  INTEGER(KIND=JPIM), PARAMETER :: iHNO3      = 5
  INTEGER(KIND=JPIM), PARAMETER :: iCH3OOH    = 6 
    ! Make copy - in TM5 chemisrty I use the other name
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O2H    = iCH3OOH 
  INTEGER(KIND=JPIM), PARAMETER :: iCH2O      = 7
  INTEGER(KIND=JPIM), PARAMETER :: iNO        = 8
  INTEGER(KIND=JPIM), PARAMETER :: iHO2       = 9
  INTEGER(KIND=JPIM), PARAMETER :: iCH3       = 10
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O      = 11
  INTEGER(KIND=JPIM), PARAMETER :: iHCO       = 12
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O2     = 13 
  INTEGER(KIND=JPIM), PARAMETER :: iOH        = 14
  INTEGER(KIND=JPIM), PARAMETER :: iNO2       = 15
  INTEGER(KIND=JPIM), PARAMETER :: iN2O5      = 16
  INTEGER(KIND=JPIM), PARAMETER :: iHO2NO2    = 17
  INTEGER(KIND=JPIM), PARAMETER :: iNO3       = 18
  INTEGER(KIND=JPIM), PARAMETER :: iN2O       = 19
  INTEGER(KIND=JPIM), PARAMETER :: iH2O       = 20
  INTEGER(KIND=JPIM), PARAMETER :: iOCLO      = 21
  INTEGER(KIND=JPIM), PARAMETER :: iHCL       = 22
  INTEGER(KIND=JPIM), PARAMETER :: iCLONO2    = 23
  INTEGER(KIND=JPIM), PARAMETER :: iHOCL      = 24
  INTEGER(KIND=JPIM), PARAMETER :: iCL2       = 25
  INTEGER(KIND=JPIM), PARAMETER :: iHBR       = 26
  INTEGER(KIND=JPIM), PARAMETER :: iBRONO2    = 27 
  INTEGER(KIND=JPIM), PARAMETER :: iCL2O2     = 28
  INTEGER(KIND=JPIM), PARAMETER :: iHOBR      = 29
  INTEGER(KIND=JPIM), PARAMETER :: iBRCL      = 30
  INTEGER(KIND=JPIM), PARAMETER :: iCFC11     = 31
  INTEGER(KIND=JPIM), PARAMETER :: iCFC12     = 32
  INTEGER(KIND=JPIM), PARAMETER :: iCFC113    = 33
  INTEGER(KIND=JPIM), PARAMETER :: iCFC114    = 34
  INTEGER(KIND=JPIM), PARAMETER :: iCFC115    = 35
  INTEGER(KIND=JPIM), PARAMETER :: iCCL4      = 36
  INTEGER(KIND=JPIM), PARAMETER :: iCLNO2     = 37
  INTEGER(KIND=JPIM), PARAMETER :: iCH3CCL3   = 38
  INTEGER(KIND=JPIM), PARAMETER :: iCH3CL     = 39
  INTEGER(KIND=JPIM), PARAMETER :: iHCFC22    = 40
  INTEGER(KIND=JPIM), PARAMETER :: iCH3BR     = 41
  INTEGER(KIND=JPIM), PARAMETER :: iHF        = 42
  INTEGER(KIND=JPIM), PARAMETER :: iHA1301    = 43
  INTEGER(KIND=JPIM), PARAMETER :: iHA1211    = 44
  INTEGER(KIND=JPIM), PARAMETER :: iCHBR3     = 45
  INTEGER(KIND=JPIM), PARAMETER :: iCLOO      = 46
  INTEGER(KIND=JPIM), PARAMETER :: iO         = 47
  INTEGER(KIND=JPIM), PARAMETER :: iO1D       = 48
  INTEGER(KIND=JPIM), PARAMETER :: iN         = 49
  INTEGER(KIND=JPIM), PARAMETER :: iCLO       = 50
  INTEGER(KIND=JPIM), PARAMETER :: iCL        = 51
  INTEGER(KIND=JPIM), PARAMETER :: iBR        = 52
  INTEGER(KIND=JPIM), PARAMETER :: iBRO       = 53
  INTEGER(KIND=JPIM), PARAMETER :: iH         = 54
  INTEGER(KIND=JPIM), PARAMETER :: iH2        = 55
  INTEGER(KIND=JPIM), PARAMETER :: iCO2       = 56
  INTEGER(KIND=JPIM), PARAMETER :: iBR2       = 57
  INTEGER(KIND=JPIM), PARAMETER :: iCH2BR2    = 58
  INTEGER(KIND=JPIM), PARAMETER :: iStratAer  = 59
  INTEGER(KIND=JPIM), PARAMETER :: iPAR       = 60      
  INTEGER(KIND=JPIM), PARAMETER :: iETH       = 61      
  INTEGER(KIND=JPIM), PARAMETER :: iOLE       = 62      
  INTEGER(KIND=JPIM), PARAMETER :: iALD2      = 63      
  INTEGER(KIND=JPIM), PARAMETER :: iPAN       = 64      
  INTEGER(KIND=JPIM), PARAMETER :: iROOH      = 65      
  INTEGER(KIND=JPIM), PARAMETER :: iORGNTR    = 66      
  INTEGER(KIND=JPIM), PARAMETER :: iISOP      = 67      
  INTEGER(KIND=JPIM), PARAMETER :: iSO2       = 68      
  INTEGER(KIND=JPIM), PARAMETER :: iDMS       = 69      
  INTEGER(KIND=JPIM), PARAMETER :: iNH3       = 70      
  INTEGER(KIND=JPIM), PARAMETER :: iSO4       = 71      
  INTEGER(KIND=JPIM), PARAMETER :: iNH4       = 72      
  INTEGER(KIND=JPIM), PARAMETER :: iMSA       = 73      
  INTEGER(KIND=JPIM), PARAMETER :: imgly      = 74      
  INTEGER(KIND=JPIM), PARAMETER :: iO3S       = 75      
  INTEGER(KIND=JPIM), PARAMETER :: iRn222     = 76      
  INTEGER(KIND=JPIM), PARAMETER :: iPb210     = 77      
  INTEGER(KIND=JPIM), PARAMETER :: iC2O3      = 78      
  INTEGER(KIND=JPIM), PARAMETER :: iROR       = 79      
  INTEGER(KIND=JPIM), PARAMETER :: iRXPAR     = 80      
  INTEGER(KIND=JPIM), PARAMETER :: iXO2       = 81      
  INTEGER(KIND=JPIM), PARAMETER :: iXO2N      = 82      
  INTEGER(KIND=JPIM), PARAMETER :: iNH2       = 83 
  INTEGER(KIND=JPIM), PARAMETER :: iPSC       = 84      
  INTEGER(KIND=JPIM), PARAMETER :: iCH3OH     = 85      
  INTEGER(KIND=JPIM), PARAMETER :: iHCOOH   =   86      
  INTEGER(KIND=JPIM), PARAMETER :: iMCOOH   =   87      
  INTEGER(KIND=JPIM), PARAMETER :: iC2H6    =   88      
  INTEGER(KIND=JPIM), PARAMETER :: iETHOH   =   89      
  INTEGER(KIND=JPIM), PARAMETER :: iC3H8    =   90      
  INTEGER(KIND=JPIM), PARAMETER :: iC3H6    =   91      
  INTEGER(KIND=JPIM), PARAMETER :: iterp    =   92      
  INTEGER(KIND=JPIM), PARAMETER :: iISPD    =   93      
  INTEGER(KIND=JPIM), PARAMETER :: iNO3_A   =   94      
  INTEGER(KIND=JPIM), PARAMETER :: iacet    =   95      
  INTEGER(KIND=JPIM), PARAMETER :: iACO2    =   96      
  INTEGER(KIND=JPIM), PARAMETER :: iIC3H7O2 =   97      
  INTEGER(KIND=JPIM), PARAMETER :: iHYPROPO2=   98      
  INTEGER(KIND=JPIM), PARAMETER :: iNOXA    =   99     
  INTEGER(KIND=JPIM), PARAMETER :: IOCS     =  100
  INTEGER(KIND=JPIM), PARAMETER :: ISO3     =  101
  INTEGER(KIND=JPIM), PARAMETER :: IH2SO4   =  102
  INTEGER(KIND=JPIM), PARAMETER :: ICH3O2NO2=  103 
  INTEGER(KIND=JPIM), PARAMETER :: IHONO    =  104 
  INTEGER(KIND=JPIM), PARAMETER :: IHCN     =  105
  INTEGER(KIND=JPIM), PARAMETER :: ICH3CN   =  106
  INTEGER(KIND=JPIM), PARAMETER :: IXYL     =  107
  INTEGER(KIND=JPIM), PARAMETER :: ITOL     =  108
  INTEGER(KIND=JPIM), PARAMETER :: IAROO2   =  109
  ! Isoprene degredation products 
  INTEGER(KIND=JPIM), PARAMETER :: IHPALD1  =  110
  INTEGER(KIND=JPIM), PARAMETER :: IHPALD2  =  111 
  INTEGER(KIND=JPIM), PARAMETER :: IISOPOOH =  112
  INTEGER(KIND=JPIM), PARAMETER :: IGLY     =  113         
  INTEGER(KIND=JPIM), PARAMETER :: IGLYALD  =  114
  INTEGER(KIND=JPIM), PARAMETER :: IHYAC    =  115
  ! Isoprene peroxy from LIM0    
  INTEGER(KIND=JPIM), PARAMETER :: IISOPBO2 =  116
  INTEGER(KIND=JPIM), PARAMETER :: IISOPDO2 =  117  

! SOA..
  INTEGER(KIND=JPIM), PARAMETER :: ISOG1    =  118      
  INTEGER(KIND=JPIM), PARAMETER :: ISOG2A   =  119      
  INTEGER(KIND=JPIM), PARAMETER :: ISOG2B   =  120    

  ! Dummy index to track jo2 budget (but it doesn't work). Keep io2 < NCHEM+3
  INTEGER(KIND=JPIM), PARAMETER :: IO2 = NTRACE + 2
  INTEGER(KIND=JPIM)            :: IACID, IAIR

  INTEGER(KIND=JPIM), PARAMETER :: ntemp=155  ! = nthigh-ntlow
  REAL(KIND=JPRB),dimension(ntrace,ntemp) :: henry   ! heterogeneous removal rates

  ! Boundary conditions at surface for stratospheric species
  ! Which currently don't feature emissions (!) 
  ! Just a simple parameter, taken from BASCOE field for 1 April 2008
  ! An exception is CHBR3, where loss in troposphere is implicitly assumed
  INTEGER(KIND=JPIM), PARAMETER  :: NBC = 19

  INTEGER(KIND=JPIM),DIMENSION(NBC), PARAMETER ::BASCOE_BC=(/ &
  & IN2O   , ICFC11 , ICFC12 , ICFC113, ICFC114, &
  & ICFC115, ICH2BR2, ICCL4 , ICH3CCL3, IHCFC22, &
  & IHA1301, IHA1211, ICH3Br, ICHBR3  , ICH3CL , &
  & ICO2   , IBRO   , IHBR,   IOCS    /)
  CHARACTER (LEN = 12), PARAMETER, DIMENSION(NBC) :: BASCOE_BCNAME = &
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
     &        "bro         " , &
     &        "hbr         " , &
     &        "ocs         " /)

  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL=(/ &
  &  3.27E-7 , 2.33E-10, 5.20E-10, 7.27E-11, 1.63E-11, &
  &  8.43E-12, 1.17E-12, 83.2E-12, 3.70E-12, 2.29E-10, &
  &  3.30E-12, 3.75E-12, 6.65E-12, 1.17E-12, 5.34E-10, &
  &  400.E-6 , 0.10E-12, 0.10E-12, 266.E-12 /)
  
  !adopt Meinshausen mean surface data (mole fraction) for 1750 / 1850
  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL_1750=(/ &
  &  2.74E-7 , 0.00E-10, 0.00E-10, 0.00E-11, 0.00E-12, &
  &  0.00E-12, 1.20E-12, 2.94E-14, 0.00E-11, 0.00E-10, &
  &  0.00E-12, 0.00E-12, 5.30E-12, 1.20E-12, 4.57E-10, &
  &  277.E-6 , 0.00E-12, 0.00E-12, 266.E-12 /)

  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL_1850=(/ &
  &  2.73E-7 , 0.00E-10, 0.00E-10, 0.00E-11, 0.00E-12, &
  &  0.00E-10, 1.20E-12, 2.94E-14, 0.00E-11, 0.00E-10, &
  &  0.00E-12, 0.00E-12, 5.30E-12, 1.20E-12, 4.57E-10, &
  &  284.E-6 , 0.00E-12, 0.00E-12, 266.E-12 /)
 

END MODULE BASCOETM5_MODULE
