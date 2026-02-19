! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CONVECTION_LAYER(YDSURF, &
 ! Input quantities
  & YDMODEL,KDIM, LDSLPHY, STATE, TENDENCY_CML, TENDENCY_DYN, PAUX, PPERT, &
 ! Input/Output quantities
  & LLKEYS, PDIAG, AUXL, PERTL, FLUX, PSURF, GEMSL, &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *CONVECTION_LAYER* - Layer routine calling full convection scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDSLPHY  : Key for SL phys
! state    : Derived variable for model state
! tendency_cml : D. V. for stored model tendencies from processes before 
! tendency_dyn : Dynamics tendencies
! PAUX     : Derived variables for general auxiliary quantities
! PPERT    : Derived variable for incoming perturbations etc... 

!     ==== Input/output ====
! LLKEYS       : Derived variable with keys
! PDIAG        : Derived variable for diagnostic quantities
! AUXL         : Derived variables for local auxiliary quantities
! PERTL        : Derived variables for local perturbartions
! FLUX         : Derived variable for fluxes
! PSURF        : Derived variables for general surface quantities
! GEMSL        : Derived variables for local GEMS quantities

!    ==== Output tendencies from convection ====
! tendency_loc :  Output process tendencies


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
!     M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!     P. Lopez (8 Dec 2020) Added separate fields for lightning parameterization.
!     M. Leutbecher (Oct 2020) SPP abstraction
!     --------------

!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, &
   & KEYS_LOCAL_TYPE, PERTURB_LOCAL_TYPE, FLUX_TYPE, GEMS_LOCAL_TYPE, &
   & PERTURB_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
LOGICAL                        , INTENT (IN)   :: LDSLPHY
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_CML
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_DYN
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (PERTURB_TYPE)            , INTENT (IN)   :: PPERT
TYPE (KEYS_LOCAL_TYPE)         , INTENT(INOUT) :: LLKEYS
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (PERTURB_LOCAL_TYPE)      , INTENT(INOUT) :: PERTL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JARP, JL
REAL(KIND=JPRB)    :: ZGP2DSPP(KDIM%KLON, YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL)  !SPP pattern

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cucalln.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CONVECTION_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSD_VN=>YDSURF%YSD_VN)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN

DO JARP=1, YDSPP_CONFIG%SM%NRFTOTAL
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZGP2DSPP(JL,JARP)=PPERT%PGP2DSPP(JL,1,JARP)
  ENDDO
ENDDO

CALL CUCALLN &
  & ( YDMODEL%YRML_PHY_RAD%YRERAD,YDMODEL%YRML_PHY_SLIN,YDMODEL%YRML_PHY_EC,YDMODEL%YRML_GCONF%YGFL, &
  & YDMODEL%YRML_CHEM%YRCHEM, YDMODEL%YRML_GCONF%YRSPP_CONFIG, &
  & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON  , KDIM%KSMAX  , KDIM%KLEV, &
  & LLKEYS%LLLAND, LDSLPHY,&
  & TSPHY,YDMODEL%YRML_PHY_G%YRVDF%RVDIFTS,&
  & STATE%T    , STATE%Q   , STATE%U   , STATE%V,   AUXL%ZLISUM,&
  & PAUX%PVERVEL, FLUX%PDIFTQ, FLUX%PDIFTS,  PAUX%PAPRS,&
  & PAUX%PRSF1  , PAUX%PRS1  , PAUX%PGEOM1, PAUX%PGEOMH, PAUX%PGAW,&
  & PERTL%ZCUCONVCA, ZGP2DSPP, &
  & TENDENCY_CML, TENDENCY_LOC, TENDENCY_DYN,& 
  & PSURF%PSD_VN(:,YSD_VN%YACPR%MP),&
  & AUXL%ITOPC  , AUXL%IBASC , PDIAG%ITYPE, &
  & PDIAG%ICBOT , PDIAG%ICTOP, AUXL%IBOTSC, LLKEYS%LLCUM , LLKEYS%LLSC,&
  & PDIAG%ICBOT_LIG , PDIAG%ICTOP_LIG, LLKEYS%LLCUM_LIG,&
  & LLKEYS%LLSHCV,&
  & AUXL%ZLCRIT_AER,&
  & PDIAG%ZLU   , PDIAG%ZLUDE, PDIAG%ZLUDELI, PDIAG%ZSNDE, PDIAG%PMFU , PDIAG%PMFD,&
  & FLUX%PDIFCQ , FLUX%PDIFCS, FLUX%PFHPCL, FLUX%PFHPCN,&
  & FLUX%PFPLCL , FLUX%PFPLCN, PDIAG%ZLRAIN,PDIAG%PRSUD,& 
  & FLUX%PSTRCU, FLUX%PSTRCV, FLUX%PFCCQL, FLUX%PFCCQN,&
  & PDIAG%PMFUDE_RATE ,    PDIAG%PMFDDE_RATE ,   PDIAG%PCAPE,  PERTL%ZWMEAN, PDIAG%PVDISCU, PERTL%ZDISSCU,&
  & GEMSL%ITRAC  , GEMSL%ZCEN  , GEMSL%ZTENC, GEMSL%ZSCAV )


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CONVECTION_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_LAYER
