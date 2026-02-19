! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE VARBC_TABLE
!-----------------------------------------------------------------------
!   VARBC_table - data structures and other definitions used across
!                 the varBC routines. But data storage is never global
!                 and belongs in VARBC_CLASS.

!   Author.
!   -------
!   Alan Geer 11-Jun-2009 - Needed to avoid circular dependencies between 
!                           varbc_pred and varbc_setup. 

!   Modifications.
!   --------------
!   Paul Poli 03-Mar-2011 - Added jpmxnqcparms max number of VARBC QC params.
!   Alan Geer 27-Jul-2015 - Objectify VarBC for OOPS
!   N. Bormann 1-Aug-2016 - CVarBC
!
!------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: JPMXNQCPARMS = 2 ! max number of VARBC QC params

! Cold start behaviour
INTEGER(KIND=JPIM), PARAMETER :: JPSTART_ZERO=0, JPSTART_PRESCRIBED=1, JPSTART_MODE=2

! Table for administration of bias groups
! ---------------------------------------
TYPE TYPE_VARBC
  INTEGER(KIND=JPIM)              :: PREVDATE     ! previous date when group was encountered
  INTEGER(KIND=JPIM)              :: NPARAM       ! number of bias parameters
  INTEGER(KIND=JPIM)              :: NCOUNT       ! data count 
  CHARACTER(LEN=8)                :: OBSCLASS     ! observation class
  CHARACTER(LEN=80)               :: GROUPKEY     ! class-dependent group description
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPREDCS(:)   ! list of predictors
  REAL(KIND=JPRB),    ALLOCATABLE :: APARAMS(:)   ! bias parameters - latest estimate
  REAL(KIND=JPRB),    ALLOCATABLE :: ZPARAMS(:)   ! bias parameters - prescribed values
  INTEGER(KIND=JPIM)              :: NCSTART      ! coldstart option
  LOGICAL                         :: LLMASKRS     ! flag for radiosonde masking
  LOGICAL                         :: LLMASKCLD    ! flag for cloud-cover masking
  LOGICAL                         :: LLMODE       ! flag for mode correction
  LOGICAL                         :: LLINCR       ! flag for incremental solution
  LOGICAL,            ALLOCATABLE :: LLCONST(:)   ! flag for constant parameters
  REAL(KIND=JPRB)                 :: OBSERR       ! rms obs error
! Parameters for CVARBC
  LOGICAL                         :: LCVARBC      ! flag for constrained VarBC
  REAL(KIND=JPRB)                 :: BIASCOR_CON  ! prescribed bias value for CVarBC
  REAL(KIND=JPRB)                 :: BIASERR_CON  ! error for the size of the bias correction
  REAL(KIND=JPRB)                 :: ALPHA_CON    ! tuning factor for CVARBC cost function
! End parameters for CVARBC
  REAL(KIND=JPRB),    ALLOCATABLE :: APARAM5(:)   ! trajectory parameter valuess
  REAL(KIND=JPRB),    ALLOCATABLE :: AGRAD(:)     ! gradient w/r to bias parameters
  REAL(KIND=JPRB),    ALLOCATABLE :: AGRAD0(:)    ! initial gradient w/r to bias parameters
  INTEGER(KIND=JPIM), ALLOCATABLE :: NHSTFGDEP(:) ! histogram of background departures
  REAL(KIND=JPRB)                 :: DFGDEP       ! histogram range
  REAL(KIND=JPRB)                 :: QCPARMS(JPMXNQCPARMS)   ! QC parameters
  INTEGER(KIND=JPIM), ALLOCATABLE :: MPREDXCNT(:,:) ! number of observations per predictor covariance
  REAL(KIND=JPRB),    ALLOCATABLE :: APREDMEAN(:)   ! mean predictor
  REAL(KIND=JPRB),    ALLOCATABLE :: APREDXCOV(:,:) ! predictor covariance
  INTEGER(KIND=JPIM)              :: NCOMP          ! local number of predictor contributions
  INTEGER(KIND=JPIM)              :: MCOMP          ! local current predictor contributions
  REAL(KIND=JPRB),    ALLOCATABLE :: APREDMEAN_COMP(:,:)   ! mean predictor contributions
  REAL(KIND=JPRB),    ALLOCATABLE :: APREDXCOV_COMP(:,:,:) ! predictor covariance contributions
  REAL(KIND=JPRB),    ALLOCATABLE :: APARAMS_COMP(:,:) ! aparams partial sum 
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPARAMS_COMP(:)   ! aparams number of terms

END TYPE TYPE_VARBC

! Table for administration of bias groups background constraint
! -------------------------------------------------------------

TYPE TYPE_VARBC_BGC
  REAL(KIND=JPRB),    ALLOCATABLE :: BKGERR(:)    ! background error std dev
  REAL(KIND=JPRB),    ALLOCATABLE :: APARAM0(:)   ! background parameter values
  REAL(KIND=JPRB),    ALLOCATABLE :: ACHGVAR(:,:) ! change-of-variable operator
  REAL(KIND=JPRB),    ALLOCATABLE :: ACVARIN(:,:) ! inverse change-of-variable
  LOGICAL                         :: LLBKGCON     ! flag for background constraint
END TYPE TYPE_VARBC_BGC

END MODULE VARBC_TABLE
