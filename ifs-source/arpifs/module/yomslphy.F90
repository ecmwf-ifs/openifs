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

MODULE YOMSLPHY

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE


TYPE :: TSLPHY

!     ------------------------------------------------------------------

! * Variables to keep inter-timestep memory in physics (split physics or VD tendencies for SC).

! RSLWX   : level of implicitness of semi-Lagrangian/physics.
! NVTEND  : third dimension of SAVTEND (number of 3D fields).
! NVTEND_VD  : third dimension of SAVTEND_VD (number of 3D fields).
! SAVTEND : buffer to store physical tendencies.
! SAVTEND5: trajectory for SAVTEND.
! SAVTEND_VD : buffer to store tendencies from vertical diffusion and radiation.

REAL(KIND=JPRB),ALLOCATABLE :: RSLWX(:)
INTEGER(KIND=JPIM) :: NVTEND, NVTEND_VD
REAL(KIND=JPRB),ALLOCATABLE :: SAVTEND(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SAVTEND5(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SAVTEND_VD(:,:,:,:)
! Pointers for SAVTEND and SAVTEND_VD
INTEGER(KIND=JPIM) :: MU_SAVTEND,MV_SAVTEND,MT_SAVTEND
! Pointers used exclusively in SAVTEND
INTEGER(KIND=JPIM) :: MUS_SAVTEND,MVS_SAVTEND,MTS_SAVTEND,MQS_SAVTEND,MSPPTCLEAR_SAVTEND
! Pointer used exclusively in SAVTEND_VD
INTEGER(KIND=JPIM) :: MQ_SAVTEND,M_RAD_SW,M_RAD_LW,M_CONV_MFD,M_CONV_MFU,M_CONV_LUDELI(4)

!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TSLPHY
!============================================================================
!TYPE(TSLPHY), POINTER :: YRSLPHY => NULL()

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TSLPHY), INTENT(IN) :: SELF
  INTEGER      , INTENT(IN) :: KDEPTH
  INTEGER      , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMSLPHY:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_g%yrslphy : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVTEND = ', SELF%NVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVTEND_VD = ', SELF%NVTEND_VD
  IF (ALLOCATED(SELF%SAVTEND)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SAVTEND allocated of shape ', &
    & SHAPE(SELF%SAVTEND)
  IF (ALLOCATED(SELF%SAVTEND5)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SAVTEND5 allocated of shape ',&
    & SHAPE(SELF%SAVTEND5)
  IF (ALLOCATED(SELF%SAVTEND_VD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SAVTEND_VD allocated of shape ', &
    & SHAPE(SELF%SAVTEND_VD)
  IF (ALLOCATED(SELF%RSLWX)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSLWX allocated of shape ',&
    & SHAPE(SELF%RSLWX)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MU_SAVTEND = ', SELF%MU_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MV_SAVTEND = ', SELF%MV_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MT_SAVTEND = ', SELF%MT_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MUS_SAVTEND = ', SELF%MUS_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MVS_SAVTEND = ', SELF%MVS_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MTS_SAVTEND = ', SELF%MTS_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MQS_SAVTEND = ', SELF%MQS_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSPPTCLEAR_SAVTEND = ', SELF%MSPPTCLEAR_SAVTEND
  IF (ALLOCATED(SELF%SAVTEND_VD)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MQ_SAVTEND = ', SELF%MQ_SAVTEND
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_RAD_SW = ', SELF%M_RAD_SW
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_RAD_LW = ', SELF%M_RAD_LW
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_CONV_MFD = ', SELF%M_CONV_MFD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_CONV_MFD = ', SELF%M_CONV_MFD
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_CONV_LUDELI = ', SELF%M_CONV_LUDELI(:)
  ENDIF
IF (LHOOK) CALL DR_HOOK('YOMSLPHY:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMSLPHY
