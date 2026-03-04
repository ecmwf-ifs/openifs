! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain constants used in UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      MODULE UKCA_CONSTANTS

      USE PARKIND1  ,ONLY : JPIM     ,JPRB

      IMPLICIT NONE
      SAVE

!!#include "c_mdi.h"
!!#include "c_sulchm.h"
!!#include "c_rmol.h"


! UKCA_ MODE constants

! MODE constants already defined in the UM
! 1) Original values
      REAL(KIND=JPRB), PARAMETER :: PPI=3.14159265358979323846_JPRB
      REAL(KIND=JPRB), PARAMETER :: AVC=6.022E23_JPRB
      REAL(KIND=JPRB), PARAMETER :: ZBOLTZ=1.3807E-23_JPRB
      REAL(KIND=JPRB), PARAMETER :: VKARMN=0.4_JPRB
      REAL(KIND=JPRB), PARAMETER :: RA=287.05_JPRB
      REAL(KIND=JPRB), PARAMETER :: GG=9.80665_JPRB
      REAL(KIND=JPRB), PARAMETER :: RAD_E=6371229.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: MM_DA=AVC*ZBOLTZ/RA
      REAL(KIND=JPRB), PARAMETER :: RHOSUL=1800.0E0_JPRB      ! UM is slightly different
!      REAL(KIND=JPRB), PARAMETER :: MMSUL=0.098E0_JPRB
      REAL(KIND=JPRB), PARAMETER :: RR=8.314_JPRB
!      PARAMETER(BCONST=3.44E13)        ! Volume_mode
!      PARAMETER(NU_H2SO4=3.0)          ! Volume_mode
!      PARAMETER(CONVERT=1.0e-18)       ! Volume_mode
      REAL(KIND=JPRB),PARAMETER :: RHOW=1000.0_JPRB ! Volume_mode
      REAL(KIND=JPRB),PARAMETER :: MMW=0.018_JPRB ! Volume_mode
!      PARAMETER(RHOSUL=1800.0)`        ! Volume_mode

! 2) UM definitions of MODE constants
!       REAL(KIND=JPRB), PARAMETER :: PPI=pi
!       REAL(KIND=JPRB), PARAMETER :: AVC=Avogadro
!       REAL(KIND=JPRB), PARAMETER :: ZBOLTZ=Boltzmann
!       REAL(KIND=JPRB), PARAMETER :: VKARMN=vKman
!       REAL(KIND=JPRB), PARAMETER :: RA=R
!       REAL(KIND=JPRB), PARAMETER :: RR=rmol
!       REAL(KIND=JPRB), PARAMETER :: GG=g
!       REAL(KIND=JPRB), PARAMETER :: RAD_E=Earth_Radius
!       REAL(KIND=JPRB), PARAMETER :: RHOSUL=rho_so4       ! 1769 or 1800 ?
!       REAL(KIND=JPRB), PARAMETER :: RHOW=rho_water       ! Water density (kg/m^3)

! 3) MODE/aerosol chemistry constants not already in UM
!      REAL(KIND=JPRB), PARAMETER :: MM_DA=AVC*ZBOLTZ/RA !
      REAL(KIND=JPRB), PARAMETER :: NMOL=1.0E2_JPRB          !
      REAL(KIND=JPRB), PARAMETER :: TDAYS=0.0_JPRB           !
      REAL(KIND=JPRB), PARAMETER :: EMS_EPS=1.0E-40_JPRB     !
      REAL(KIND=JPRB), PARAMETER :: CONC_EPS=1.0E-8_JPRB     !
      REAL(KIND=JPRB), PARAMETER :: DN_EPS=1.0E-40_JPRB      !
      REAL(KIND=JPRB), PARAMETER :: BCONST=3.44E13_JPRB      ! Volume_mode
      REAL(KIND=JPRB), PARAMETER :: NU_H2SO4=3.0_JPRB        ! Volume_mode
      REAL(KIND=JPRB), PARAMETER :: CONVERT=1.0E-18_JPRB     ! Volume_mode
      REAL(KIND=JPRB), PARAMETER :: H_PLUS=1.0E-5_JPRB       ! cloud/rain H+ concentration

! molecular masses in kg/mol

      LOGICAL, PARAMETER  ::  L_UKCA_DIURNAL_ISOPEMS   = .TRUE.
      REAL(KIND=JPRB), PARAMETER :: M_AIR        = 0.02897_JPRB       ! Air
      REAL(KIND=JPRB), PARAMETER :: MMSUL        = 0.09808_JPRB       ! H2SO4
!      REAL(KIND=JPRB), PARAMETER :: mmw          = 0.0180154     ! H2O


!  -------------------------------------------------------------------
!  Conversion factor from vmr to mmr for each species.
!             vmr*c_species = mmr
!             c_species = m_species/m_air  (m_air = 28.97)
!
!  Followed by molecular masses for each species (e.g. m_ch4) in g/mol
!  -------------------------------------------------------------------

      REAL(KIND=JPRB), PARAMETER :: C_O3P        = 0.5523_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_O1D        = 0.5523_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_O3         = 1.657_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_NO         = 1.036_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_NO3        = 2.140_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_NO2        = 1.588_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_N2O5       = 3.728_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HO2NO2     = 2.727_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HONO2      = 2.175_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HNO3       = C_HONO2   ! used in nitrate.F90
      REAL(KIND=JPRB), PARAMETER :: C_OH         = 0.5868_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HO2        = 1.139_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H2         = 0.06904_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H2O2       = 1.174_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CH4        = 0.5523_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_C          = 0.4142_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CO         = 0.9665_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CO2        = 1.5188_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HCHO       = 1.036_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEOO       = 1.622_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H2O        = 0.6213_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H2OS       = 0.6213_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEOOH      = 1.657_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HONO       = 1.622_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_O2         = 1.105_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_N2         = 0.9665_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_C2H6       = 1.036_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ETOO       = 2.106_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ETOOH      = 2.140_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECHO      = 1.519_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_TOTH       = 1.000_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECO3      = 2.589_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_PAN        = 4.177_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_C3H8       = 1.519_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_PROO       = 2.589_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_PROOH      = 2.623_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ETCHO      = 2.002_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ETCO3      = 3.072_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ME2CO      = 2.002_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECOCH2OO  = 3.072_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECOCH2OOH = 3.107_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_PPAN       = 4.660_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEONO2     = 2.658_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_N          = 0.48325_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H          = 0.03452_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_N2O        = 1.5188_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CFCL3      = 4.7480_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2CL2     = 4.1783_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CLO        = 1.7784_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HCL        = 1.2604_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CLONO2     = 3.3668_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HOCL       = 1.8129_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_OCLO       = 2.3309_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_BRO        = 3.315_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_BRONO2     = 4.9034_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HBR        = 2.7970_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HOBR       = 3.3495_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_BRCL       = 3.9884_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEBR       = 3.2805_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_SO2        = 2.2112_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_SO3        = 2.7615_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ME2S       = 2.145_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_DMS        = 2.145_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_DMSO       = 2.6965_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_OCS        = 2.0711_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_COS        = 2.0711_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H2S        = 1.1766_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CS2        = 2.6282_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_SAD        = 4.1255_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MSA        = 3.317_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_S          = 1.1046_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_H2SO4      = 3.385_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2CLCFCL2 = 6.4722_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CHF2CL     = 2.9858_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECCL3     = 4.6082_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CCL4       = 5.3158_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECL       = 1.7432_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2CLBR    = 5.7128_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF3BR      = 5.1432_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CL         = 1.2261_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CL2O2      = 3.5568_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_BR         = 2.7627_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CH2BR2     = 6.0013_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECF2CL    = 3.4673_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2BR2     = 7.2489_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2BRCF2BR = 8.9748_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2CLCF3   = 5.3314_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_CF2CLCF2CL = 5.8992_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECFCL2    = 4.0352_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_C5H8       = 2.3473_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ISO2       = 4.0387_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ISOOH      = 4.0732_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ISON       = 5.3504_JPRB  ! Might be revised to 5.0052
!                       for RAQ chem where ISON = (NO3)C4H6CHO: C5H7NO4: 145
      REAL(KIND=JPRB), PARAMETER :: C_MACR       = 2.4163_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MACRO2     = 4.1077_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MACROOH    = 4.1422_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MPAN       = 5.0742_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HACET      = 2.5544_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MGLY       = 2.4853_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_NALD       = 3.6244_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HCOOH      = 1.5878_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECO3H     = 2.6234_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MECO2H     = 2.0711_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_NH3        = 0.5879_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MONOTERP   = 4.7034_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_SEC_ORG    = 5.1782_JPRB    ! Molecular weight=150.


!     Extra species for RAQ chemistry
      REAL(KIND=JPRB), PARAMETER :: C_C4H10      = 2.0021_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_C2H4       = 0.9665_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_C3H6       = 1.4498_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_RNC2H4     = 3.6244_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_RNC3H6     = 4.1077_JPRB
!     rnc2h4 & rnc3h6 are CH2(NO3)CHO & CH3CH(NO3)CHO
      REAL(KIND=JPRB), PARAMETER :: C_CH3OH      = 1.1046_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEOH       = 1.1046_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_TOLUENE    = 3.1757_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_OXYLENE    = 3.6590_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEMALD     = 3.3828_JPRB
!     MEMALDIAL is CH3-CO-CH=CH-CHO, ring fragmentation_JPRB
!     product from degradation of toluene and o-xylene
      REAL(KIND=JPRB), PARAMETER :: C_BUOOH      = 3.1067_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEK        = 2.4853_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MVK        = 2.4163_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MVKOOH     = 4.1422_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_GLY        = 2.0021_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_RNC5H8     = 5.0052_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ORGNIT     = 5.5230_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_BUOO       = 3.0721_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEKO2      = 3.5554_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HOC2H4O2   = 2.6579_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HOC3H6O2   = 3.1412_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HOIPO2     = 4.0387_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_TOLP1      = 3.7280_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_OXYL1      = 4.2113_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEMALD1    = 3.9351_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_HOMVKO2    = 4.1077_JPRB

!     Extra species for EXTTC chemistry
      REAL(KIND=JPRB), PARAMETER :: C_ISOO       = 4.0387_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MACROO     = 4.1077_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_APIN       = 4.9645_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_PROPEOO    = 2.5198_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_PROPEOOH   = 2.5544_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ONITU      = 3.5554_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEKOO      = 3.5554_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_MEKOOH     = 3.5889_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ETEOO      = 2.0366_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ALKA       = 2.0021_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ALKAOO     = 3.0721_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ALKAOOH    = 3.1067_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_AROM       = 3.4173_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_AROMOO     = 4.4874_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_AROMOOH    = 4.5219_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_BSVOC1     = 4.9645_JPRB   ! as APIN
      REAL(KIND=JPRB), PARAMETER :: C_BSVOC2     = 4.9645_JPRB   ! as APIN
      REAL(KIND=JPRB), PARAMETER :: C_B2NDRY     = 4.9645_JPRB   ! as APIN
      REAL(KIND=JPRB), PARAMETER :: C_ASVOC1     = 3.4173_JPRB   ! as AROM
      REAL(KIND=JPRB), PARAMETER :: C_ASVOC2     = 3.4173_JPRB   ! as AROM
      REAL(KIND=JPRB), PARAMETER :: C_A2NDRY     = 3.4173_JPRB   ! as AROM
      REAL(KIND=JPRB), PARAMETER :: C_BSOA       = 5.1778_JPRB   ! 150.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ASOA       = 5.1778_JPRB   ! 150.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: C_ISOSVOC1   = 2.3473_JPRB   ! as C5H8
      REAL(KIND=JPRB), PARAMETER :: C_ISOSVOC2   = 2.3473_JPRB   ! as C5H8
      REAL(KIND=JPRB), PARAMETER :: C_ISOSOA     = 4.4874_JPRB   ! 130.0_JPRB

!     molecular masses in g/mol of emitted species,
!     for budget calculations

      REAL(KIND=JPRB), PARAMETER :: M_HO2     =  33.007_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CH4     =  16._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CO      =  28._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_HCHO    =  30._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_C2H6    =  30._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_C3H8    =  44._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MECHO   =  44._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_NO2     =  46._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_N2O5    = 108.01_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ME2CO   =  58._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ISOP    =  68._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_NO      =  30._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_C       =  12._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MONOTERP =  136.24_JPRB

!     molecular masses of stratospheric species, for which surface
!     mmrs are prescribed

      REAL(KIND=JPRB), PARAMETER :: M_HCL        =  36.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_N2O        =  44._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CLO        =  51.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_HOCL       =  52.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_OCLO       =  67.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CLONO2     =  97.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CF2CL2     = 121._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CFCL3      = 137.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_HBR        =  81._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MEBR       =  95._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BRO        =  96._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_HOBR       =  97._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BRCL       = 115.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BRONO2     = 142._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CF2CLCFCL2 = 187.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CHF2CL     =  86.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MECCL3     = 133.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CCL4       = 154._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MECL       =  50.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CF2CLBR    = 165.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CF3BR      = 149._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CH2BR2     = 173.835_JPRB

! sulphur containing, etc.
      REAL(KIND=JPRB), PARAMETER :: M_OCS        =  60._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_COS        =  60._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_H2S        =  34.086_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CS2        =  76.14_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_DMS        =  62.1_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_DMSO       =  78.13_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ME2S       =  62.1_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MSA        =  96.1_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_SEC_ORG    =  150.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_S          =  32.07_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_SO2        =  64.06_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_SO3        =  80.06_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_SO4        =  96.06_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_H2SO4      =  98.07_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_NH3        =  17.03_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_NH42SO4    =  132.16_JPRB

      REAL(KIND=JPRB), PARAMETER :: M_CL         = 35.5_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CL2O2      = 103._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BR         = 80._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_H2         = 2.016_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_H2O        = 18.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MECOCH2OOH = 90.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ISOOH      = 118.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MPAN       = 147.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_PPAN       = 135.0_JPRB

!     Extra masses for RAQ or other chemistries
      REAL(KIND=JPRB), PARAMETER :: M_C5H8    =  68.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_C4H10   =  58.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_C2H4    =  28.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_C3H6    =  42.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_TOLUENE =  92.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_OXYLENE = 106.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_CH3OH   =  32.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MEOH    =  32.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BUOOH   =  90.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MVKOOH  = 120.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ORGNIT  = 160.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MACROOH = 120.0_JPRB

      REAL(KIND=JPRB), PARAMETER :: M_HONO = 47.0_JPRB
!     Extra masses for Wesely scheme
      REAL(KIND=JPRB), PARAMETER :: M_MACR    = 70.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ETCHO   = 58.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_NALD    = 105.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MGLY    = 72.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_HACET  = 74.0_JPRB


!     Extra masses for EXTTC chemistry
      REAL(KIND=JPRB), PARAMETER :: M_APIN     =  136._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MVK      =  70._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MEK      =  72._JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ALKA     =  58._JPRB        ! as butane
      REAL(KIND=JPRB), PARAMETER :: M_AROM     =  99._JPRB        ! (toluene + xylene)/2
      REAL(KIND=JPRB), PARAMETER :: M_BSVOC1   = 144.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BSVOC2   = 144.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ASVOC1   = 99.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ASVOC2   = 99.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ISOSVOC1 = 68.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ISOSVOC2 = 68.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ONITU    = 102.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_BSOA     = 150.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ASOA     = 150.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ISOSOA   = 130.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_ALKAOOH  = 90.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_AROMOOH  = 130.0_JPRB
      REAL(KIND=JPRB), PARAMETER :: M_MEKOOH   = 104.0_JPRB

!     The mass of organic nitrate is an approximation,
!     calculated as the average of ORGNIT formed by two
!     reacs. in UKCA_CHEMCO_RAQ:
!      NO2 + TOLP1 --> ORGNIT (A)
!      NO2 + OXYL1 --> ORGNIT (B)
!      * TOL  = methylbenzene       = C6H5(CH3)
!        OXYL = 1,2-dimethylbenzene = C6H4(CH3)2
!      * TOL  + OH --> TOLP1: C6H4(OH)(CH3)  = methyl phenol
!        OXYL + OH --> OXYL1: C6H3(OH)(CH3)2 = dimethyl phenol
!      * ORGNIT A: TOLP1 + NO2 ~ C6H3(CH3)(OH)NO2  ~
!                  C7H7NO3: methyl nitrophenol   -> 153
!        ORGNIT B: OXYL1 + NO2 ~ C6H2(CH3)2(OH)NO2 ~
!                  C8H9NO3: dimethyl nitrophenol -> 167
!  -------------------------------------------------------------------

      END MODULE UKCA_CONSTANTS
