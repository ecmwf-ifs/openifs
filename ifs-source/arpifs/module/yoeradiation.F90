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

MODULE YOERADIATION
  
USE PARKIND1,         ONLY :   JPRB,JPIM
USE radiation_config, ONLY :   config_type

IMPLICIT NONE

SAVE


! This derived type contains configuration information for the
! radiation scheme plus a few additional variables and parameters
! needed for the IFS interface to it
TYPE :: TRADIATION

  ! Configuration information for the radiation scheme
  type(config_type)  :: rad_config

  ! Ultraviolet weightings
  INTEGER(KIND=JPIM) :: NWEIGHT_UV
  INTEGER(KIND=JPIM) :: IBAND_UV(100)
  REAL(KIND=JPRB)    :: WEIGHT_UV(100)
  ! Photosynthetically active radiation weightings
  INTEGER(KIND=JPIM) :: NWEIGHT_PAR
  INTEGER(KIND=JPIM) :: IBAND_PAR(100)
  REAL(KIND=JPRB)    :: WEIGHT_PAR(100)
  ! Background aerosol is specified in an ugly way: using the old
  ! Tegen fields that are in terms of optical depth, and converted to
  ! mass mixing ratio via the relevant mass-extinction
  ! coefficient. The following are the indices to the aerosol types
  ! used to describe tropospheric and stratospheric background
  ! aerosol, assigned later.
  INTEGER(KIND=JPIM) :: ITYPE_TROP_BG_AER  = 0  ! Hydrophobic organic
  INTEGER(KIND=JPIM) :: ITYPE_STRAT_BG_AER = 0 ! Stratospheric sulphate
  ! Mass-extinction coefficient (m2 kg-1) of tropospheric and
  ! stratospheric background aerosol at 550 nm
  REAL(KIND=JPRB)    :: TROP_BG_AER_MASS_EXT  = 0.0_JPRB
  REAL(KIND=JPRB)    :: STRAT_BG_AER_MASS_EXT = 0.0_JPRB
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TRADIATION

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TRADIATION),  INTENT(IN) :: SELF
  INTEGER(KIND=JPIM), INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM), INTENT(IN) :: KOUTNO
  
  INTEGER(KIND=JPIM) :: IDEPTHLOC
  
  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yradiation : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // '** we should print content of rad_config, not done yet**'
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NWEIGHT_UV = ', SELF%NWEIGHT_UV
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'IBAND_UV SUM = ', SUM(SELF%IBAND_UV)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WEIGHT_UV SUM = ', SUM(SELF%WEIGHT_UV)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NWEIGHT_PAR = ', SELF%NWEIGHT_PAR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'IBAND_PAR SUM = ', SUM(SELF%IBAND_PAR)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'WEIGHT_PAR SUM = ', SUM(SELF%WEIGHT_PAR)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TROP_BG_AER_MASS_EXT = ', SELF%TROP_BG_AER_MASS_EXT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'STRAT_BG_AER_MASS_EXT = ', SELF%STRAT_BG_AER_MASS_EXT
  
END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOERADIATION
