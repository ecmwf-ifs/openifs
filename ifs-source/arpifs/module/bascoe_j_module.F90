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

MODULE BASCOE_J_MODULE

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!**   DESCRIPTION 
!     ----------
!
!   Module defining strato photolysis rates.
!
!     MODIFICATIONS.
!     --------------
!     2017-12-07: YVES CHRISTOPHE (YC) *BIRA*
!       extracted from bascoe_module.F90
!----------------

! Photolysis reaction rates (J) (chem sb14a,sb14b,sb15b)
!
  INTEGER(KIND=JPIM), PARAMETER :: ndiss = 55   ! Nb of J files (<= nb J processes)

  INTEGER(KIND=JPIM), PARAMETER :: J_O2_O       =  1
  INTEGER(KIND=JPIM), PARAMETER :: J_O3_O       =  2
  INTEGER(KIND=JPIM), PARAMETER :: J_O3_O1D     =  3
  INTEGER(KIND=JPIM), PARAMETER :: J_HO2        =  4
  INTEGER(KIND=JPIM), PARAMETER :: J_H2O2_OH    =  5
  INTEGER(KIND=JPIM), PARAMETER :: J_NO2        =  6
  INTEGER(KIND=JPIM), PARAMETER :: J_NO3_O      =  7
  INTEGER(KIND=JPIM), PARAMETER :: J_NO3_O2     =  8
  INTEGER(KIND=JPIM), PARAMETER :: J_N2O5       =  9
  INTEGER(KIND=JPIM), PARAMETER :: J_HNO3       = 10
  INTEGER(KIND=JPIM), PARAMETER :: J_HO2NO2_HO2 = 11
  INTEGER(KIND=JPIM), PARAMETER :: J_Cl2        = 12
  INTEGER(KIND=JPIM), PARAMETER :: J_Br2        = 13
  INTEGER(KIND=JPIM), PARAMETER :: J_OClO       = 14
  INTEGER(KIND=JPIM), PARAMETER :: J_Cl2O2      = 15
  INTEGER(KIND=JPIM), PARAMETER :: J_HOCl       = 16
  INTEGER(KIND=JPIM), PARAMETER :: J_ClONO2_Cl  = 17
  INTEGER(KIND=JPIM), PARAMETER :: J_BrCl       = 18
  INTEGER(KIND=JPIM), PARAMETER :: J_BrONO2_Br  = 19
  INTEGER(KIND=JPIM), PARAMETER :: J_CH2O_HCO   = 20
  INTEGER(KIND=JPIM), PARAMETER :: J_CH2O_CO    = 21
  INTEGER(KIND=JPIM), PARAMETER :: J_CH3OOH     = 22
  INTEGER(KIND=JPIM), PARAMETER :: J_HO2NO2_OH  = 23
  INTEGER(KIND=JPIM), PARAMETER :: J_ClONO2_ClO = 24
  INTEGER(KIND=JPIM), PARAMETER :: J_BrO        = 25
  INTEGER(KIND=JPIM), PARAMETER :: J_CCl4       = 26
  INTEGER(KIND=JPIM), PARAMETER :: J_CFC11      = 27
  INTEGER(KIND=JPIM), PARAMETER :: J_CFC113     = 28
  INTEGER(KIND=JPIM), PARAMETER :: J_CFC114     = 29
  INTEGER(KIND=JPIM), PARAMETER :: J_CFC115     = 30
  INTEGER(KIND=JPIM), PARAMETER :: J_CFC12      = 31
  INTEGER(KIND=JPIM), PARAMETER :: J_CH3Br      = 32
  INTEGER(KIND=JPIM), PARAMETER :: J_CH3CCl3    = 33
  INTEGER(KIND=JPIM), PARAMETER :: J_CH3Cl      = 34
  INTEGER(KIND=JPIM), PARAMETER :: J_CH4        = 35
  INTEGER(KIND=JPIM), PARAMETER :: J_CHBr3      = 36
  INTEGER(KIND=JPIM), PARAMETER :: J_ClNO2      = 37
  INTEGER(KIND=JPIM), PARAMETER :: J_ClOO       = 38
  INTEGER(KIND=JPIM), PARAMETER :: J_ClO        = 39
  INTEGER(KIND=JPIM), PARAMETER :: J_H2O        = 40
  INTEGER(KIND=JPIM), PARAMETER :: J_HA1211     = 41
  INTEGER(KIND=JPIM), PARAMETER :: J_HA1301     = 42
  INTEGER(KIND=JPIM), PARAMETER :: J_HCFC22     = 43
  INTEGER(KIND=JPIM), PARAMETER :: J_HCl        = 44
  INTEGER(KIND=JPIM), PARAMETER :: J_HOBr       = 45
  INTEGER(KIND=JPIM), PARAMETER :: J_N2O        = 46
  INTEGER(KIND=JPIM), PARAMETER :: J_NO         = 47
  INTEGER(KIND=JPIM), PARAMETER :: J_O2_O1D     = 48
  INTEGER(KIND=JPIM), PARAMETER :: J_BrONO2_BrO = 49
  INTEGER(KIND=JPIM), PARAMETER :: J_CO2        = 50
  INTEGER(KIND=JPIM), PARAMETER :: J_CH2Br2     = 51
  INTEGER(KIND=JPIM), PARAMETER :: J_H2O2_HO2   = 52
  INTEGER(KIND=JPIM), PARAMETER :: J_OCS        = 53
  INTEGER(KIND=JPIM), PARAMETER :: J_H2SO4      = 54
  INTEGER(KIND=JPIM), PARAMETER :: J_SO3        = 55

  ! VH - change capitals to lower case
  CHARACTER (LEN = 12), PARAMETER, DIMENSION(ndiss) :: jnames = &
     &      (/"o2-o        " , &
     &        "o3-o        " , &
     &        "o3-o1d      " , &
     &        "ho2         " , &
     &        "h2o2-oh     " , &
     &        "no2         " , &
     &        "no3-o       " , &
     &        "no3-o2      " , &
     &        "n2o5        " , &
     &        "hno3        " , &
     &        "ho2no2-ho2  " , &
     &        "cl2         " , &
     &        "br2         " , &
     &        "oclo        " , &
     &        "cl2o2       " , &
     &        "hocl        " , &
     &        "clono2-cl   " , &
     &        "brcl        " , &
     &        "brono2-br   " , &
     &        "ch2o-hco    " , &
     &        "ch2o-co     " , &
     &        "ch3ooh      " , &
     &        "ho2no2-oh   " , &
     &        "clono2-clo  " , &
     &        "bro         " , &
     &        "ccl4        " , &
     &        "cfc11       " , &
     &        "cfc113      " , &
     &        "cfc114      " , &
     &        "cfc115      " , &
     &        "cfc12       " , &
     &        "ch3br       " , &
     &        "ch3ccl3     " , &
     &        "ch3cl       " , &
     &        "ch4         " , &
     &        "chbr3       " , &
     &        "clno2       " , &
     &        "cloo        " , &
     &        "clo         " , &
     &        "h2o         " , &
     &        "ha1211      " , &
     &        "ha1301      " , &
     &        "hcfc22      " , &
     &        "hcl         " , &
     &        "hobr        " , &
     &        "n2o         " , &
     &        "no          " , &
     &        "o2-o1d      " , &
     &        "brono2-bro  " , &    
     &        "co2         " , &
     &        "ch2br2      " , &
     &        "h2o2-ho2    " , &
     &        "ocs         " , &
     &        "h2so4       " , &
     &        "so3         "  /)


END MODULE BASCOE_J_MODULE
