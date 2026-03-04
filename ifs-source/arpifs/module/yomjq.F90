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

MODULE YOMJQ

!   Purpose
!   -------
!     Data structures to store the components of model error
!     covariance matrix.

!   Author
!   ------
!     Yannick Tremolet

!   Modifications
!   -------------
!     Original    16-Jan-04
!     M. Chrust   3-Jan-2020 Define ERRMOD_STRUCT/OOPS cleaning
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPRB, JPIM
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD

IMPLICIT NONE
SAVE

LOGICAL         :: LJQZERO   ! Model error without a Jq term (dangerous)
LOGICAL         :: LSTATMERR ! Generate data for model error statistics
REAL(KIND=JPRB) :: FJQCOST   ! Model error term of cost function

! Model error covariance config
TYPE ERRMOD_CONFIG_STRUCT
  LOGICAL                         :: LJQHCOR            = .FALSE.      ! Use horizontal correlations
  LOGICAL                         :: LJQVCOR            = .FALSE.      ! Use vertical   correlations
  LOGICAL                         :: LSCALERR           = .FALSE.      ! Scale model error
  INTEGER(KIND=JPIM)              :: NLEVERR0, NLEVERR1                ! Limit model error to stratosphere
  INTEGER(KIND=JPIM)              :: NDIM_MODERR        = 0            !
  INTEGER(KIND=JPIM)              :: NSERR              = 0            ! Spectral truncation for model error
  INTEGER(KIND=JPIM)              :: NVE3D              = 0            ! Number of model error 3D fields
  INTEGER(KIND=JPIM)              :: NVE2D              = 0            ! Number of model error 2D fields
  INTEGER(KIND=JPIM), ALLOCATABLE :: NGRBMERRCTL(:)                    ! Grib code model error
  REAL(KIND=JPRB)                 :: TSTEP              = 0.0_JPRB
  REAL(KIND=JPRB)                 :: TSTEP_ERR          = 0.0_JPRB
  REAL(KIND=JPRB)                 :: ALPHAQ             = 0.0_JPRB     ! Coef. model error term in cost function
  REAL(KIND=JPRB)                 :: ALPHAQ_T           = 0.0_JPRB     ! Coef. temperature model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_D           = 0.0_JPRB     ! Coef. divergence  model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_V           = 0.0_JPRB     ! Coef. vorticity   model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_Q           = 0.0_JPRB     ! Coef. humidity    model error term
  REAL(KIND=JPRB)                 :: ALPHAQ_O           = 0.0_JPRB     ! Coef. ozone       model error term
END TYPE ERRMOD_CONFIG_STRUCT

! ERRMODMN: Mean profiles of model error
TYPE ERRMODMN_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE :: V2D(:)
  REAL(KIND=JPRB), ALLOCATABLE :: V3D(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: V2DINV(:)
  REAL(KIND=JPRB), ALLOCATABLE :: V3DINV(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: VALPHA(:)
  REAL(KIND=JPRB) :: VSURF
END TYPE ERRMODMN_STRUCT

! ERRHCOR: Horizontal model error correlation matrix (inverse)
TYPE ERRMODHCOR_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE :: SP2D(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: SP3D(:,:,:)
END TYPE ERRMODHCOR_STRUCT

! ERRVCOR: Vertical model error correlation matrix
TYPE ERRMODVCOR_STRUCT
  REAL(KIND=JPRB), ALLOCATABLE :: VOR(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: DIV(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: PT(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: Q(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: O3(:,:,:)
END TYPE ERRMODVCOR_STRUCT

! Model error covariance matrix
TYPE ERRMOD_STRUCT
  INTEGER(KIND=JPIM), ALLOCATABLE :: NVALERR(:)         ! Info local Ns
  TYPE (ERRMOD_CONFIG_STRUCT)     :: CONFIG             ! Configuration
  TYPE (ERRMODMN_STRUCT)          :: ERRMODMN           ! Mean profiles of model error
  TYPE (ERRMODHCOR_STRUCT)        :: ERRHCOR            ! Horizontal model error correlation matrix (inverse)
  TYPE (ERRMODVCOR_STRUCT)        :: ERRVCOR, ERRVCORIN ! Vertical model error correlation matrix
END TYPE ERRMOD_STRUCT

TYPE(ERRMOD_STRUCT) :: YGERRMOD

END MODULE YOMJQ
