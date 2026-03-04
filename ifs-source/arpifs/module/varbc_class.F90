! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE VARBC_CLASS
!-----------------------------------------------------------------------
!   VARBC_table - Variational bias correction class

!   Author.
!   -------
!   Alan Geer 27-Jul-2015 - Represent VarBC as an object, for OOPS

!                         - Design note: much of the old VarBC still
!                           lives below as independent modules (when
!                           this process is finished, there will be no
!                           global data in them) but varbc should only
!                           be accessed through this class

!   Modifications.
!   --------------
!   Niels Bormann  1-Aug-2016  CVarBC
!   V. Guidard  06-Oct-2016  Adapt for MERGE_VARBC program
!
!------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB, JPRD
USE YOMCT0, ONLY : L_OOPS
USE VARBC_TABLE, ONLY : TYPE_VARBC,TYPE_VARBC_BGC
USE VARBC_PRED, ONLY : JPREDNAME

USE GOM_PLUS , ONLY : TYPE_GOM_PLUS
USE OML_MOD  , ONLY : OML_MAX_THREADS
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMANCS , ONLY : RMDI
USE YOMLUN, ONLY : NULOUT
USE DBASE_MOD

IMPLICIT NONE
PRIVATE

TYPE, PUBLIC :: CLASS_VARBC
  PRIVATE

  INTEGER(KIND=JPIM), PUBLIC :: NOBGROUP   ! number of bias groups
  INTEGER(KIND=JPIM), PUBLIC :: NOBPARAM   ! total number of bias parameters (dimension of augmented control vector)
  INTEGER(KIND=JPIM), PUBLIC :: MXNPRED    ! number of predictors
  INTEGER(KIND=JPIM)         :: MXNPARAM   ! number of params (just MXNPRED+1)
  INTEGER(KIND=JPIM), PUBLIC :: MXBODYPRED ! number of predictors with body-dependence
  INTEGER(KIND=JPIM)         :: BODYPREDSTART ! start ID of predictors with body-dependence
  REAL(KIND=JPRB),    PUBLIC :: AFJPCOST = 0.0_JPRB  ! VarBC cost function; initialised to initial contribution

  TYPE(TYPE_VARBC), ALLOCATABLE, PUBLIC :: YVARBC(:) ! main VarBC table
  TYPE(TYPE_VARBC_BGC), ALLOCATABLE, PUBLIC :: YVARBC_BGC(:)
  CHARACTER(LEN=JPREDNAME), ALLOCATABLE :: CPREDDESC(:) ! predictor descriptions
  REAL(KIND=JPRB), ALLOCATABLE :: PPREDNORM(:,:) ! predictor normalization: global mean, stdv
  INTEGER(KIND=JPIM), ALLOCATABLE :: IHSTFGDEP_COMP(:,:,:) ! nhstfgdep partial sum
  LOGICAL, ALLOCATABLE :: L_SET_HAS_CVARBC(:) ! For OOPS CVARBC optimization (indicates obs_set has obs with LCVARBC=True)
END TYPE

CONTAINS

!------------------------------------------------------------------------
!------------------------------------------------------------------------
END MODULE VARBC_CLASS
