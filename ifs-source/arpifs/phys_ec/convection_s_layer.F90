! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CONVECTION_S_LAYER(YDSURF, YDMODEL, &
 ! Input quantities
  & KDIM, STATE, TENDENCY_CML, &
  & PTENDENCY_VD9, PAUX,&
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
! state    : Derived variable for  model state
! tendency_cml : D. V. for model resulting tendencies
! PTENDENCY_VD9 : tendencies from previous timestep vertical diffusion
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
!     F. Vana       (Jan 2023): new seq. physics order

!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX   , ONLY : TSURF
USE TYPE_MODEL           , ONLY : MODEL
USE PARKIND1             , ONLY : JPRB
USE YOMHOOK              , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   &                  AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, &
   &                  KEYS_LOCAL_TYPE, FLUX_TYPE, GEMS_LOCAL_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)               , INTENT(INOUT) :: YDSURF
TYPE(MODEL)               , INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)     , INTENT (IN)   :: KDIM
TYPE (STATE_TYPE)         , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)         , INTENT (IN)   :: TENDENCY_CML
REAL(KIND=JPRB), TARGET   , INTENT (IN)   :: PTENDENCY_VD9(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND_VD)
TYPE (AUX_TYPE)           , INTENT (IN)   :: PAUX
TYPE (KEYS_LOCAL_TYPE)    , INTENT(INOUT) :: LLKEYS
TYPE (AUX_DIAG_TYPE)      , INTENT(INOUT) :: PDIAG
TYPE (AUX_DIAG_LOCAL_TYPE), INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)          , INTENT(INOUT) :: FLUX
TYPE (SURF_AND_MORE_TYPE) , INTENT(INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)    , INTENT(INOUT) :: GEMSL
TYPE (STATE_TYPE)         , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------

REAL(KIND=JPRB),POINTER :: ZTENDENCY_VD9(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cucalln2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CONVECTION_S_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDML_PHY_SLIN=>YDMODEL%YRML_PHY_SLIN,& 
 & YDML_PHY_EC=>YDMODEL%YRML_PHY_EC,YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY, &
 & YDVDF=>YDMODEL%YRML_PHY_G%YRVDF,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, YSD_VN=>YDSURF%YSD_VN, &
 & MT_SAVTEND=>YDSLPHY%MT_SAVTEND, MU_SAVTEND=>YDSLPHY%MU_SAVTEND, &
 & MV_SAVTEND=>YDSLPHY%MV_SAVTEND, MQ_SAVTEND=>YDSLPHY%MQ_SAVTEND, &
 & NVTEND_VD=>YDSLPHY%NVTEND_VD, RVDIFTS=>YDVDF%RVDIFTS)
!     ------------------------------------------------------------------

!  Set vertical diffusion tendencies from previous time-step
ALLOCATE(ZTENDENCY_VD9(KDIM%KLON,KDIM%KLEV,NVTEND_VD))
ZTENDENCY_VD9(:,:,:)=0._JPRB

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN2

CALL CUCALLN2 &
  & (YDMODEL%YRML_PHY_EC%YRTHF,YDMODEL%YRCST, YDERAD,YDML_PHY_SLIN,YDML_PHY_EC, &
  & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON  , KDIM%KLEV,&
  & LLKEYS%LLLAND, LLKEYS%LLSLPHY, LLKEYS%LLRAIN1D, &
  & TSPHY,RVDIFTS,&
  & STATE%T     , STATE%Q  , STATE%U    , STATE%V,&
  & PAUX%PVERVEL, FLUX%PDIFTQ(:,KDIM%KLEV+1), FLUX%PDIFTS(:,KDIM%KLEV+1), PAUX%PAPRS,&
  & PAUX%PRSF1  , PAUX%PRS1  , PAUX%PGEOM1, PAUX%PGEOMH, PAUX%PGAW,&
  & TENDENCY_LOC%T, TENDENCY_CML%T, ZTENDENCY_VD9(:,:,MT_SAVTEND), &
  & TENDENCY_LOC%Q, TENDENCY_CML%Q, ZTENDENCY_VD9(:,:,MQ_SAVTEND), &
  & TENDENCY_LOC%U, TENDENCY_CML%U, ZTENDENCY_VD9(:,:,MU_SAVTEND), &
  & TENDENCY_LOC%V ,TENDENCY_CML%V, ZTENDENCY_VD9(:,:,MV_SAVTEND), &
  & PSURF%PSD_VN(:,YSD_VN%YACPR%MP),&
  & AUXL%ITOPC  , AUXL%IBASC , PDIAG%ITYPE,&
  & PDIAG%ICBOT , PDIAG%ICTOP, AUXL%IBOTSC, LLKEYS%LLCUM , LLKEYS%LLSC,&
  & PDIAG%ZLU   , PDIAG%ZLUDE, PDIAG%PMFU , PDIAG%PMFD,&
  & FLUX%PDIFCQ , FLUX%PDIFCS, FLUX%PFHPCL, FLUX%PFHPCN,&
  & FLUX%PFPLCL , FLUX%PFPLCN, FLUX%PSTRCU, FLUX%PSTRCV, FLUX%PFCCQL, FLUX%PFCCQN,&
  & PDIAG%PMFUDE_RATE ,    PDIAG%PMFDDE_RATE ,   PDIAG%PCAPE,&
  & GEMSL%ITRAC  , GEMSL%ZCEN  , GEMSL%ZTENC,  GEMSL%ZSCAV )

! Cleaning
DEALLOCATE(ZTENDENCY_VD9)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CONVECTION_S_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_S_LAYER
