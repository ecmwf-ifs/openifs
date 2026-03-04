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
MODULE YOERDI

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERDI* - COEFFICIENTS WITHIN RADIATION INTERFACE
!     -----------------------------------------------------------------

TYPE :: TERDI
REAL(KIND=JPRB) :: RRAE
REAL(KIND=JPRB) :: RSUNDUR
REAL(KIND=JPRB) :: RCARDI
REAL(KIND=JPRB) :: RCH4
REAL(KIND=JPRB) :: RN2O
REAL(KIND=JPRB) :: RNO2
REAL(KIND=JPRB) :: RO3
REAL(KIND=JPRB) :: RCCL4
REAL(KIND=JPRB) :: RCFC11
REAL(KIND=JPRB) :: RCFC12
REAL(KIND=JPRB) :: RCFC22
REAL(KIND=JPRB) :: REPCLC
REAL(KIND=JPRB) :: REPH2O
REAL(KIND=JPRB) :: RCCO2, RCCH4, RCN2O, RCNO2, RCCFC11, RCCFC12, RCCFC22, RCCCL4
REAL(KIND=JPRB) :: RSOLINC
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TERDI
!============================================================================

!!TYPE(TERDI), POINTER :: YRERDI => NULL()

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     Original  J.-J. MORCRETTE       E.C.M.W.F.      89/07/14
!     Modified  P. Viterbo    99/03/26    Surface tiling
!     Modified  P. Viterbo    24/05/2004  surf library
!     Modified JJMorcrette    2005/01/19  GHG and Solar constant variability

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! RRAE   : EFFECT OF EARTH'S CURVATURE ON COSINE SOLAR ZENITH ANGLE
! RSUNDUR: MINIMUM DIRECT SOLAR FOR COMPUTING SOLAR DURATION
!
! RCARDI, RCH4, RN2O, RNO2, RO3, RCFC11, RCFC12, RCFC22, RCCL4
!    Mass mixing ratios of trace gases (CARDI=Carbon Dioxide), updated
!    each timestep in UPDRGAS (if LHGHG=TRUE).  These can be thought
!    of as annual-mean concentrations at the surface, and are used to
!    scale the monthly mean latitude-pressure climatologies
!
! RCCO2, RCCH4, RCN2O, RCCFC11, RCCFC12, RCCFC22, RCCCL4
!    If LHGHG=FALSE then instead we assume the gas concentrations are
!    constant with time and are taken from these *volume* mixing
!    ratios, which are set in SUECRAD and may be overridden in the
!    NAERAD namelist.
!
! REPCLC : SECURITY TO AVOID ZERO OR ONE CLOUD COVERS
! REPH2O : SECURITY TO AVOID WATER VAPOUR CONTENT IN A LAYER
!          TO BE MORE THAN THE RESPECTIVE VALUE AT SATURATION.
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

IMPLICIT NONE
CLASS(TERDI), INTENT(IN) :: SELF
INTEGER     , INTENT(IN) :: KDEPTH
INTEGER     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOERDI:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yrerdi : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRAE = ', SELF%RRAE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSUNDUR = ', SELF%RSUNDUR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCARDI = ', SELF%RCARDI
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCH4 = ', SELF%RCH4
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RN2O = ', SELF%RN2O
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RNO2 = ', SELF%RNO2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RO3 = ', SELF%RO3
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCL4 = ', SELF%RCCL4
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCFC11 = ', SELF%RCFC11
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCFC12 = ', SELF%RCFC12
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCFC22 = ', SELF%RCFC22
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPCLC = ', SELF%REPCLC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPH2O = ', SELF%REPH2O
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCO2 = ', SELF%RCCO2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCH4 = ', SELF%RCCH4
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCN2O = ', SELF%RCN2O
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCNO2 = ', SELF%RCNO2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCFC11 = ', SELF%RCCFC11
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCFC12 = ', SELF%RCCFC12
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCFC22 = ', SELF%RCCFC22
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCCL4 = ', SELF%RCCCL4
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSOLINC = ', SELF%RSOLINC
IF (LHOOK) CALL DR_HOOK('YOERDI:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOERDI
