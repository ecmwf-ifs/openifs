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

MODULE YOEAERSNK

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERSNK* - CONTROL OPTIONS FOR AEROSOLS' SINKS
!     ------------------------------------------------------------------

TYPE :: TEAERSNK
INTEGER(KIND=JPIM):: NDRYDEPVEL_DYN
INTEGER(KIND=JPIM) :: NDRYDEP
REAL(KIND=JPRB) :: R_R, R_S
REAL(KIND=JPRB) :: RAERTS(21)
REAL(KIND=JPRB) :: RFRAER

REAL(KIND=JPRB) :: RRHTAB(12)
REAL(KIND=JPRB) :: RSSGROWTH_RHTAB(12)
REAL(KIND=JPRB) :: RSSDENS_RHTAB(12)
REAL(KIND=JPRB) :: RMMD_DD(9)
REAL(KIND=JPRB) :: RRHO_DD(9)
REAL(KIND=JPRB) :: RMMD_NI(2)
REAL(KIND=JPRB) :: RRHO_NI(2)
REAL(KIND=JPRB) :: RMMD_SS(9)
REAL(KIND=JPRB) :: RRHO_SS(9)
REAL(KIND=JPRB) :: RHO_WAT, RHO_ICE
REAL(KIND=JPRB) :: RHAMAKER

REAL(KIND=JPRB) :: RSO2CV1, RSO2CV2, RSUCV1 , RSUCV2 , RVSO2CV1, RVSO2CV2
!---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEAERSNK
!=====================================================================

TYPE(TEAERSNK), POINTER :: YREAERSNK => NULL()

!     ------------------------------------------------------------------
! NDRYDEP    : dry deposition 1: a la GEMS; 2: "exp" (a la D_GRG_4.6)
! R_R        : mean radius for rain drops (m)
! R_S        : mean radius for snow crystals (m)
! RAERTS     : inverse time scale (s-1) - used only in AER_LOSS{,_TL,_AD}, never called
! re-evaporation constants
! RFRAER     : for aerosols
! RRHTAB     : reference relative humidity in growth factor index for aerosol optical properties

! coefficient for parametrisation conversion SO2 to SO4
! RSO2CV1 and RSO2CV2 as used in aer_so2so4
! RSUCV1  and RSUCV2  possibly read in namelist naeaer 
!     -----------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

IMPLICIT NONE
CLASS(TEAERSNK), INTENT(IN) :: SELF
INTEGER        , INTENT(IN) :: KDEPTH
INTEGER        , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOEAERSNK:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_aer%yreaersnk : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDRYDEPVEL_DYN = ', SELF%NDRYDEPVEL_DYN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'R_R = ', SELF%R_R
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'R_S = ', SELF%R_S
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAERTS SUM = ', SUM(SELF%RAERTS)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RFRAER = ', SELF%RFRAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSSGROWTH_RHTAB SUM = ', SUM(SELF%RSSGROWTH_RHTAB)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSSDENS_RHTAB SUM = ', SUM(SELF%RSSDENS_RHTAB)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMMD_DD SUM = ', SUM(SELF%RMMD_DD)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRHO_DD SUM = ', SUM(SELF%RRHO_DD)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMMD_NI SUM = ', SUM(SELF%RMMD_NI)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRHO_NI SUM = ', SUM(SELF%RRHO_NI)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMMD_SS SUM = ', SUM(SELF%RMMD_SS)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRHO_SS SUM = ', SUM(SELF%RRHO_SS)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RHO_WAT = ', SELF%RHO_WAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RHO_ICE = ', SELF%RHO_ICE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RHAMAKER = ', SELF%RHAMAKER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSO2CV1 = ', SELF%RSO2CV1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSO2CV2 = ', SELF%RSO2CV2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSUCV1 = ', SELF%RSUCV1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSUCV2 = ', SELF%RSUCV2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVSO2CV1 = ', SELF%RVSO2CV1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVSO2CV2 = ', SELF%RVSO2CV2
IF (LHOOK) CALL DR_HOOK('YOEAERSNK:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEAERSNK

