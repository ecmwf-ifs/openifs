! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CONVECTION_S_LAYER(YDSURF, &
 ! Input quantities
  & YDERAD,YDML_PHY_SLIN,YDML_PHY_EC,YDPHY2,KDIM, LDSLPHY, LDRAIN1D, STATE, TENDENCY_CML, PAUX, PVDIFTS,&
 ! Input/Output quantities
  & LLKEYS, PDIAG, AUXL, FLUX, PSURF, GEMSL, &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *CONVECTION_LAYER* - Layer routine calling simplified convection scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDSLPHY  : Key for SL phys
! LDRAIN1D : Key for 1D rain assim
! state    : Derived variable for  model state
! tendency_cml : D. V. for model resulting tendencies
! PAUX     : Derived variables for general auxiliary quantities

!     ==== Input/output ====
! LLKEYS       : Derived variable with keys
! PDIAG        : Derived variable for diagnostic quantities
! AUXL         : Derived variables for local auxiliary quantities
! FLUX         : Derived variable for fluxes
! PSURF        : Derived variables for general surface quantities
! GEMSL        : Derived variables for local GEMS quantities

!    ==== Output tendencies from convection ====
! tendency_loc :  Derived variables with process tendencies


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-22  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     F. Vana  Oct-2013  Bugfix

!-----------------------------------------------------------------------

USE MODEL_PHYSICS_ECMWF_MOD      , ONLY : MODEL_PHYSICS_ECMWF_TYPE
USE MODEL_PHYSICS_SIMPLINEAR_MOD , ONLY : MODEL_PHYSICS_SIMPLINEAR_TYPE
USE YOERAD                       , ONLY : TERAD
USE SURFACE_FIELDS_MIX           , ONLY : TSURF
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   &                  AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, &
   &                  KEYS_LOCAL_TYPE, FLUX_TYPE, GEMS_LOCAL_TYPE
USE YOMPHY2  , ONLY : TPHY2

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TERAD), INTENT(INOUT) :: YDERAD
TYPE(MODEL_PHYSICS_ECMWF_TYPE),INTENT(INOUT)     :: YDML_PHY_EC
TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE),INTENT(INOUT):: YDML_PHY_SLIN
TYPE(TPHY2)               , INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE)     , INTENT (IN)   :: KDIM
LOGICAL                   , INTENT (IN)   :: LDSLPHY
LOGICAL                   , INTENT (IN)   :: LDRAIN1D
TYPE (STATE_TYPE)         , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)         , INTENT (IN)   :: TENDENCY_CML
TYPE (AUX_TYPE)           , INTENT (IN)   :: PAUX
REAL(KIND=JPRB)           , INTENT (IN)   :: PVDIFTS
TYPE (KEYS_LOCAL_TYPE)    , INTENT(INOUT) :: LLKEYS
TYPE (AUX_DIAG_TYPE)      , INTENT(INOUT) :: PDIAG
TYPE (AUX_DIAG_LOCAL_TYPE), INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)          , INTENT(INOUT) :: FLUX
TYPE (SURF_AND_MORE_TYPE) , INTENT(INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)    , INTENT(INOUT) :: GEMSL
TYPE (STATE_TYPE)         , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cucalln2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CONVECTION_S_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSD_VN=>YDSURF%YSD_VN)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN2

CALL CUCALLN2 &
  & ( YDERAD,YDML_PHY_SLIN,YDML_PHY_EC, &
  & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON  , KDIM%KLEV,&
  & LLKEYS%LLLAND, LDSLPHY, LDRAIN1D, &
  & TSPHY,PVDIFTS,&
  & STATE%T     , STATE%Q  , STATE%U    , STATE%V,&
  & PAUX%PVERVEL, FLUX%PDIFTQ, FLUX%PDIFTS, PAUX%PAPRS,&
  & PAUX%PRSF1  , PAUX%PRS1  , PAUX%PGEOM1, PAUX%PGEOMH, PAUX%PGAW,&
  & TENDENCY_LOC%T, TENDENCY_CML%T, &
  & TENDENCY_LOC%Q, TENDENCY_CML%Q, &
  & TENDENCY_LOC%U, TENDENCY_CML%U, &
  & TENDENCY_LOC%V ,TENDENCY_CML%V, &
  & PSURF%PSD_VN(:,YSD_VN%YACPR%MP),&
  & AUXL%ITOPC  , AUXL%IBASC , PDIAG%ITYPE,&
  & PDIAG%ICBOT , PDIAG%ICTOP, AUXL%IBOTSC, LLKEYS%LLCUM , LLKEYS%LLSC,&
  & PDIAG%ZLU   , PDIAG%ZLUDE, PDIAG%PMFU , PDIAG%PMFD,&
  & FLUX%PDIFCQ , FLUX%PDIFCS, FLUX%PFHPCL, FLUX%PFHPCN,&
  & FLUX%PFPLCL , FLUX%PFPLCN, FLUX%PSTRCU, FLUX%PSTRCV, FLUX%PFCCQL, FLUX%PFCCQN,&
  & PDIAG%PMFUDE_RATE ,    PDIAG%PMFDDE_RATE ,   PDIAG%PCAPE,&
  & GEMSL%ITRAC  , GEMSL%ZCEN  , GEMSL%ZTENC,  GEMSL%ZSCAV )


!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CONVECTION_S_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_S_LAYER
