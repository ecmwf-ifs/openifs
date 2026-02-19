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

MODULE YOMMODERRMOD

!     Purpose.
!     --------
!       Controls for model error in the model

!     Author.
!     -------
!       M. Chrust

!     Modifications.
!     --------------
!       Original    03-Jan-2020
! ------------------------------------------------------------------

USE PARKIND1,ONLY : JPIM,JPRB

IMPLICIT NONE

SAVE

!     --------------------------------------------------------------------------------
!     Model error configuration
!     --------------------------------------------------------------------------------

TYPE :: TMODERR

  LOGICAL                         :: LSCALERR      =.FALSE. ! Scale model error
  LOGICAL                         :: LTAPERERR     =.FALSE. ! Taper model error
  INTEGER(KIND=JPIM)              :: NTYPE_MODERR  =0       ! Type of model error
                                                            ! 0 - No model error
                                                            ! 1 - Full state
                                                            ! 2 - Forcing (constant or periodic)
                                                            ! 3 - Bias-like error
  INTEGER(KIND=JPIM)              :: NSTEP_MODERR  =0       ! Number of intervals  in model error
  INTEGER(KIND=JPIM)              :: NORDER_MODERR =0       ! Order of model error forcing
  INTEGER(KIND=JPIM)              :: NPRTMODERR    =0       ! Printing level for diagnostics
  INTEGER(KIND=JPIM)              :: TAPERERR_INI  =0       ! Start of tapering
  INTEGER(KIND=JPIM)              :: TAPERERR_END  =0       ! End of tapering
  REAL(KIND=JPRB)                 :: TSTEP_ERR     =0.0     ! Forcing time-step
  INTEGER(KIND=JPIM), ALLOCATABLE :: MSTEPERR(:)            ! Index model error steps
  REAL(KIND=JPRB), ALLOCATABLE    :: TSTEPERR(:)            ! Taper coefficient
  REAL(KIND=JPRB), ALLOCATABLE    :: FSTEPERR(:,:)          ! Coef. periodic forcing

CONTAINS

  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

END TYPE TMODERR

!     --------------------------------------------------------------------------------
CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TMODERR), INTENT(IN) :: SELF
  INTEGER       , INTENT(IN) :: KDEPTH
  INTEGER       , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%tmoderr : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTYPE_MODERR  = ', SELF%NTYPE_MODERR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NORDER_MODERR = ', SELF%NORDER_MODERR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPRTMODERR    = ', SELF%NPRTMODERR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LTAPERERR     = ', SELF%LTAPERERR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPERERR_INI  = ', SELF%TAPERERR_INI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPERERR_END  = ', SELF%TAPERERR_END

  WRITE(KOUTNO,*) ''

END SUBROUTINE PRINT_CONFIGURATION

!     --------------------------------------------------------------------------------

END MODULE YOMMODERRMOD
