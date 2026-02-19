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

! * Variables for split ECMWF physics.

! RSLWX   : level of implicitness of semi-Lagrangian/physics.
! NVTEND  : third dimension of SAVTEND (number of 3D fields).
! SAVTEND : buffer to store physical tendencies.
! SAVTEND5: trajectory for SAVTEND.

REAL(KIND=JPRB),ALLOCATABLE :: RSLWX(:)
INTEGER(KIND=JPIM) :: NVTEND
REAL(KIND=JPRB),ALLOCATABLE :: SAVTEND(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SAVTEND5(:,:,:,:)
! Pointers for SAVTEND
INTEGER(KIND=JPIM) :: MU_SAVTEND,MV_SAVTEND,MT_SAVTEND
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
  IF (ALLOCATED(SELF%SAVTEND)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SAVTEND allocated of shape ', &
    & SHAPE(SELF%SAVTEND)
  IF (ALLOCATED(SELF%SAVTEND5)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SAVTEND5 allocated of shape ',&
    & SHAPE(SELF%SAVTEND5)
  IF (ALLOCATED(SELF%RSLWX)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSLWX allocated of shape ',&
    & SHAPE(SELF%RSLWX)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MU_SAVTEND = ', SELF%MU_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MV_SAVTEND = ', SELF%MV_SAVTEND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MT_SAVTEND = ', SELF%MT_SAVTEND
IF (LHOOK) CALL DR_HOOK('YOMSLPHY:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMSLPHY
