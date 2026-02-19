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

MODULE YOMMODERRCONF

!     Purpose.
!     --------
!       Controls for model error state/increment in 4D-Var.

!     Author.
!     -------
!       M. Chrust

!     Modifications.
!     --------------
!       Original    03-Jan-2020
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

TYPE TMODERR_CONF

  LOGICAL :: LSCALERR    =.FALSE.                          ! Scale model error for I/O (from Model config)
  LOGICAL :: LMERRPTRSET =.FALSE.                          ! True if model error pointers are set

  INTEGER(KIND=JPIM) :: NCOMP_MODERR  = 0                  ! Number of components in model error
  INTEGER(KIND=JPIM) :: NSTEP_MODERR  = 0                  ! Number of intervals  in model error
  INTEGER(KIND=JPIM) :: NDIM_MODERR   = 0                  ! Total number of components in model error
  INTEGER(KIND=JPIM) :: NTYPE_MODERR  = 0                  ! Type of model error (from Model config)
  INTEGER(KIND=JPIM) :: NORDER_MODERR = 0                  ! Order of model error forcing
  INTEGER(KIND=JPIM) :: NPRTMODERR    = 0                  ! Diagnostics

  INTEGER(KIND=JPIM) :: NSERR  = 0                         ! Spectral truncation for model error
  INTEGER(KIND=JPIM) :: NVE3D  = 0                         ! Number of model error 3D fields
  INTEGER(KIND=JPIM) :: NVE2D  = 0                         ! Number of model error 2D fields
  INTEGER(KIND=JPIM) :: NVE3SP = 0                         ! Number of model error 3D SP fields
  INTEGER(KIND=JPIM) :: NVE2SP = 0                         ! Number of model error 2D SP fields
  INTEGER(KIND=JPIM) :: NVE3GP = 0                         ! Number of model error 3D GP fields
  INTEGER(KIND=JPIM) :: NVE2GP = 0                         ! Number of model error 2D GP fields

  INTEGER(KIND=JPIM), ALLOCATABLE :: NGRBMERRCTL(:)        ! Grib code model error
  INTEGER(KIND=JPIM), ALLOCATABLE :: NGRBMERRSP(:)         ! Grib spectral model error
  INTEGER(KIND=JPIM), ALLOCATABLE :: NGRBMERRGP(:)         ! Grib gridpoint model error

  INTEGER(KIND=JPIM), ALLOCATABLE :: MERRGFL(:)            ! Index model error in GFL
  INTEGER(KIND=JPIM), ALLOCATABLE :: MERRGMV(:)            ! Index model error in GMV
  INTEGER(KIND=JPIM), ALLOCATABLE :: MERRGMVS(:)           ! Index model error in GMVS
  INTEGER(KIND=JPIM), ALLOCATABLE :: MERROBS(:)            ! Index model error in COBS
  INTEGER(KIND=JPIM), ALLOCATABLE :: MERRSP2(:)            ! Index model error in SPA2
  INTEGER(KIND=JPIM), ALLOCATABLE :: MERRSP3(:)            ! Index model error in SPA3

  REAL(KIND=JPRB)                 :: ALPHAQ    = 0.0_JPRB  ! Coef. model error term in cost function
  REAL(KIND=JPRB)                 :: ALPHAQ_T  = 0.0_JPRB  ! Coef. temperature model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_D  = 0.0_JPRB  ! Coef. divergence  model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_V  = 0.0_JPRB  ! Coef. vorticity   model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_Q  = 0.0_JPRB  ! Coef. humidity    model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_O  = 0.0_JPRB  ! Coef. ozone       model error term

  REAL(KIND=JPRB)                 :: TSTEP_ERR = 0.0_JPRB  ! (from Model config)

CONTAINS

  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

END TYPE TMODERR_CONF

!     --------------------------------------------------------------------------------
CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KOUTNO)
  IMPLICIT NONE
  CLASS(TMODERR_CONF), INTENT(IN) :: SELF
  INTEGER             , INTENT(IN) :: KOUTNO

  WRITE(KOUTNO,*) 'TMODERR_CONF : '
  WRITE(KOUTNO,*) 'LSCALERR         = ', SELF%LSCALERR
  WRITE(KOUTNO,*) 'NCOMP_MODERR     = ', SELF%NCOMP_MODERR
  WRITE(KOUTNO,*) 'NSTEP_MODERR     = ', SELF%NSTEP_MODERR
  WRITE(KOUTNO,*) 'NDIM_MODERR      = ', SELF%NDIM_MODERR
  WRITE(KOUTNO,*) 'NTYPE_MODERR     = ', SELF%NTYPE_MODERR
  WRITE(KOUTNO,*) 'NPRTMODERR       = ', SELF%NPRTMODERR
  WRITE(KOUTNO,*) 'NSERR            = ', SELF%NSERR
  WRITE(KOUTNO,*) 'NVE3D            = ', SELF%NVE3D
  WRITE(KOUTNO,*) 'NVE2D            = ', SELF%NVE2D
  WRITE(KOUTNO,*) 'NVE3SP           = ', SELF%NVE3SP
  WRITE(KOUTNO,*) 'NVE2SP           = ', SELF%NVE2SP
  WRITE(KOUTNO,*) 'NVE3GP           = ', SELF%NVE3GP
  WRITE(KOUTNO,*) 'NVE2GP           = ', SELF%NVE2GP
  WRITE(KOUTNO,*) 'TSTEP_ERR        = ', SELF%TSTEP_ERR
  WRITE(KOUTNO,*) 'ASLPHAQ          = ', SELF%ALPHAQ
  WRITE(KOUTNO,*) 'ASLPHAQ_T        = ', SELF%ALPHAQ_T
  WRITE(KOUTNO,*) 'ASLPHAQ_D        = ', SELF%ALPHAQ_D
  WRITE(KOUTNO,*) 'ASLPHAQ_V        = ', SELF%ALPHAQ_V
  WRITE(KOUTNO,*) 'ASLPHAQ_Q        = ', SELF%ALPHAQ_Q
  WRITE(KOUTNO,*) 'ASLPHAQ_O        = ', SELF%ALPHAQ_O

  WRITE(KOUTNO,*) ''

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMMODERRCONF
