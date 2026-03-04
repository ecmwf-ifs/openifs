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

MODULE YOMCOM

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!-----------------------------------------------------------------------
!**** *SUMCC* * - MODULE INCLUDING VARIABLES FOR THE COM COUPLING

!     AUTHOR.
!     -------
!        Jean-Philippe Piedelievre CNRM/GMGEC/EAC

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 99-02-24
!       03-04-08  Philippe Piedelievre CNRM/GMGEC/EAC
!-----------------------------------------------------------------------

TYPE :: TCOM

! NVCOM     : number of parameters for the COM scheme.

INTEGER(KIND=JPIM) :: NVCOM

!*    DIAGNOSTICS POUR COUPLAGE COM

! OMLDTH    PROFONDEUR DE LA COUCHE DE MELANGE OCEANIQUE
! GTTLIN    GRADIENT THERMIQUE DE LA THERMOCLINE
! SSTPRE    SST DYNAMIQUE
! SSTMSK    MASQUE DE MIXAGE DYNAMIQUE-STATISTIQUE
! LOMLDTH   ??? (missing comment)

REAL(KIND=JPRB),ALLOCATABLE:: OMLDTH(:)
REAL(KIND=JPRB),ALLOCATABLE:: GTTLIN(:)
REAL(KIND=JPRB),ALLOCATABLE:: SSTPRE(:)
REAL(KIND=JPRB),ALLOCATABLE:: SSTMSK(:)

LOGICAL :: LOMLDTH

! TRAFLX         : correction flux for transports.
REAL(KIND=JPRB),ALLOCATABLE:: TRAFLX(:)
 !---------------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
 
END TYPE TCOM
!======================================================================
!!TYPE(TCOM), POINTER :: YRCOM => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TCOM), INTENT(IN) :: SELF
INTEGER(KIND=JPIM), INTENT(IN) :: KDEPTH
INTEGER(KIND=JPIM), INTENT(IN) :: KOUTNO

INTEGER(KIND=JPIM) :: IDEPTHLOC 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMCOM:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
IDEPTHLOC = KDEPTH + 2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_aoc%yrcom : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVCOM = ', SELF%NVCOM
IF (ALLOCATED(SELF%OMLDTH)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'OMLDTH allocated of shape ',SHAPE(SELF%OMLDTH), &
 &       ' sum ',SUM(SELF%OMLDTH)
IF (ALLOCATED(SELF%GTTLIN)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GTTLIN allocated of shape ',SHAPE(SELF%GTTLIN), &
 &       ' sum ',SUM(SELF%GTTLIN)
IF (ALLOCATED(SELF%SSTPRE)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SSTPRE allocated of shape ',SHAPE(SELF%SSTPRE), &
 &       ' sum ',SUM(SELF%SSTPRE)
IF (ALLOCATED(SELF%SSTMSK)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SSTMSK allocated of shape ',SHAPE(SELF%SSTMSK), &
 &       ' sum ',SUM(SELF%SSTMSK)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LOMLDTH = ', SELF%LOMLDTH
IF (ALLOCATED(SELF%TRAFLX)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRAFLX allocated of shape ',SHAPE(SELF%TRAFLX), &
 &       ' sum ',SUM(SELF%TRAFLX)

IF (LHOOK) CALL DR_HOOK('YOMCOM:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)
END SUBROUTINE PRINT_CONFIGURATION

!     -----------------------------------------------------------------
END MODULE YOMCOM
