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

MODULE YOMCUMFS
USE PARKIND1, ONLY : JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!       ----------------------------------------------------------------
!**    ** *YOMCUMFS* CONTAINS VARIABLES FOR SIMPLIFIED CONVECTION SCHEME 
!       ----------------------------------------------------------------

TYPE :: TCUMFS
LOGICAL :: LECUMFS
LOGICAL :: LREGCV
LOGICAL :: LMFCFL2_SHSTEP
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TCUMFS
!============================================================================

!!TYPE(TCUMFS), POINTER :: YRCUMFS => NULL()

!*    *YOMCUMFS* CONTAINS VARIABLES FOR SIMPLIFIED CONVECTION SCHEME 

!     P. LOPEZ      E.C.M.W.F.    09/01/2002

!     NAME            TYPE     DESCRIPTION
!     ----            ----     -----------

!***  LECUMFS         LOGICAL  SWITCH ON SIMPLIFIED CONVECTIVE MASS-FLUX SCHEME IN TRAJECTORY
!***  LREGCV          LOGICAL  SWITCH ON REGULARIZATIONS OF SIMPLIFIED CONVECTION SCHEME 
!***                           IN TANGENT-LINEAR AND ADJOINT CALCULATIONS
!***  LMFCFL2_SHSTEP  LOGICAL  SWITCH ON REDUCTION OF MASS FLUX CFL CRITERION FOR SHORT TIME STEPS 

!     ------------------------------------------------------------------
CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

  IMPLICIT NONE
  CLASS(TCUMFS), INTENT(IN) :: SELF
  INTEGER      , INTENT(IN) :: KDEPTH
  INTEGER      , INTENT(IN) :: KOUTNO
  
  INTEGER :: IDEPTHLOC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOMCUMFS:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
  IDEPTHLOC = KDEPTH+2

  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LECUMFS = ', SELF%LECUMFS
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LREGCV = ', SELF%LREGCV
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LMFCFL2_SHSTEP = ', SELF%LMFCFL2_SHSTEP
  IF (LHOOK) CALL DR_HOOK('YOMCUMFS:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)
  
END SUBROUTINE PRINT_CONFIGURATION
END MODULE YOMCUMFS
