! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BACKSCATTER_LAYER(&
 ! Input quantities
 &  YDDIM, YDML_PHY_STOCH,YDDYN,YDPHY2,KDIM, PAUX, STATE_T0, STATE_TMP,&
 ! Input/Output quantities
 &  PPERT, PERTL, PDIAG, PSURF, TENDENCY)

!**** *BACKSCATTER_LAYER* - Layer routine calling the backscatter scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities
! state_t0 : Derived variable for model state at t0
! state_tmp: Derived variable for model state entering the convection scheme

!     ==== Input/output ====
! PPERT    : Derived variable for in/out perturbartions
! PERTL    : Derived variable for local perturbartions
! PDIAG    : Derived variable for diagnostics quantities
! PSURF    : Derived variable for surface (and more) fields
! tendency : local tendency


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!       Note: Output directly applied to the tendencies which are then IN/OUT here

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-23  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE MODEL_PHYSICS_STOCHAST_MOD , ONLY : MODEL_PHYSICS_STOCHAST_TYPE
USE YOMDYN   , ONLY : TDYN
USE YOMPHY2  , ONLY : TPHY2
USE YOMDIM   , ONLY : TDIM
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE,&
  &  PERTURB_TYPE, PERTURB_LOCAL_TYPE, SURF_AND_MORE_TYPE,&
  &  AUX_DIAG_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM)                     , INTENT (IN)   :: YDDIM
TYPE(TDYN)                      ,INTENT(INOUT) :: YDDYN
TYPE(MODEL_PHYSICS_STOCHAST_TYPE),INTENT(INOUT):: YDML_PHY_STOCH
TYPE(TPHY2)                     ,INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE_T0
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE_TMP
TYPE (PERTURB_TYPE)            , INTENT(INOUT) :: PPERT
TYPE (PERTURB_LOCAL_TYPE)      , INTENT(INOUT) :: PERTL
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY


!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "spbsgpupd.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('BACKSCATTER_LAYER',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN

CALL SPBSGPUPD(YDDIM,YDML_PHY_STOCH,YDDYN,YDPHY2,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,PPERT%PVORTGRADX,PPERT%PVORTGRADY,&
         & PPERT%PVORT,TENDENCY%U,TENDENCY%V,PAUX%PAPRSF,STATE_T0%T,PDIAG%PMFUDE_RATE,PDIAG%PMFU,PERTL%ZWMEAN,&
         & STATE_TMP%T,STATE_TMP%Q,PAUX%PRSF1,PPERT%PSTREAM,PPERT%PTEMP,PPERT%PSTOPHCA,&
         & PPERT%PTOTDISS,PPERT%PTOTDISS_SMOOTH,PERTL%ZDISSGW,PERTL%ZDISSCU,PSURF%PSD_XA,PSURF%PSD_X2,&
         & KDIM%KLEVX,KDIM%KFLDX,KDIM%KFLDX2,PAUX%PGEMU,STATE_T0%U,STATE_T0%V,PAUX%PAPRS,PDIAG%ITYPE)


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BACKSCATTER_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE BACKSCATTER_LAYER
