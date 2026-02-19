! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE CLOUD_LAYER( &
 ! Input quantities
  & YDSURF, YDECLDP,YDECUMF,YDEPHLI,YDPHY2,YDERAD,YDEPHY,YDVDF,YDSPP_CONFIG,YGFL,&
  & KDIM, LDSLPHY, PAUX, PPERT, STATE, & 
  & TENDENCY_CML, TENDENCY_DYN, TENDENCY_VDF, PRAD, &
  & PSURF, LLKEYS, &
 ! Input/Output quantities
  & AUXL, FLUX, PDIAG, FSD, &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *CLOUD_LAYER* - Layer routine calling cloud scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDSLPHY  : key to activate SL physics
! PAUX     : Derived variables for general auxiliary quantities
! PPERT    : Derived variable for incoming perturbations etc... 
! state    : Derived variable for model state
! tendency_cml : D. V. for model tendencies  from processes before 
! tendency_dyn : D. V. for model tendencies from explicit dynamics
! tendency_vdf : D. V. for model tendencies from turbulence scheme
! PRAD     : D. V. for radiative quantities
! PSURF    : D.V. for surface quantities
! LLKEYS   : D.V. for local switches

!     ==== Input/output ====
! AUXL         : Derived variable for local quantites
! FLUX         : Derived variable for fluxes
! PDIAG        : Derived variable for diagnostics
! PFSD         : cloud heterogeneity variable
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
!      Original : 2012-11-28  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      2015-01-10  R. Forbes  Precip type diags now here so called every timestp
!      M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!      M Ahlgrimm 2017-11-11 add cloud heterogeneity FSD
!      R. Forbes  2020-11-15 Remove TENDENCY_TMP and add various inputs to cloudsc
!      M. Leutbecher 2020-10-12 SPP abstraction
!-----------------------------------------------------------------------

USE YOECLDP            , ONLY : TECLDP, NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE YOECUMF            , ONLY : TECUMF
USE YOEPHLI            , ONLY : TEPHLI
USE YOERAD             , ONLY : TERAD
USE YOEPHY             , ONLY : TEPHY
USE YOEVDF             , ONLY : TVDF
USE SPP_MOD            , ONLY : TSPP_CONFIG
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE VARIABLE_MODULE    , ONLY : VARIABLE_3D
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMPHYDER          , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, PERTURB_TYPE, &
 &                              FLUX_TYPE, AUX_DIAG_LOCAL_TYPE, AUX_RAD_TYPE, KEYS_LOCAL_TYPE, &
 &                              SURF_AND_MORE_TYPE, AUX_DIAG_TYPE
USE YOM_YGFL           , ONLY : TYPE_GFLD
USE YOMPHY2            , ONLY : TPHY2

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TECLDP),INTENT(INOUT) :: YDECLDP
TYPE(TECUMF),INTENT(INOUT) :: YDECUMF
TYPE(TEPHLI),INTENT(INOUT) :: YDEPHLI
TYPE(TPHY2) ,INTENT(INOUT) :: YDPHY2
TYPE(TERAD) ,INTENT(INOUT) :: YDERAD
TYPE(TEPHY) ,INTENT(INOUT) :: YDEPHY
TYPE(TVDF)  ,INTENT(INOUT) :: YDVDF
TYPE(TSPP_CONFIG) ,INTENT(IN) :: YDSPP_CONFIG
TYPE(TYPE_GFLD),INTENT(INOUT) :: YGFL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
LOGICAL                        , INTENT (IN)   :: LDSLPHY
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (PERTURB_TYPE)            , INTENT (IN)   :: PPERT
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_CML
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_DYN
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_VDF
TYPE (AUX_RAD_TYPE)            , INTENT (IN)   :: PRAD
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (KEYS_LOCAL_TYPE)         , INTENT (IN)   :: LLKEYS
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (VARIABLE_3D)             , INTENT(INOUT) :: FSD
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JRF, JL, JK
REAL(KIND=JPRB)    :: ZGP2DSPP(KDIM%KLON, YDSPP_CONFIG%SM%NRFTOTAL)  !SPP pattern
REAL(KIND=JPRB)    :: ZFSD(KDIM%KLON,KDIM%KLEV)


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cloudsc.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUD_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSD_VF=>YDSURF%YSD_VF, YSD_VD=>YDSURF%YSD_VD)

! ------------------------------------------------------------------------------
!
!*       1. Initialisation
!
! ------------------------------------------------------------------------------
DO JRF=1, YDSPP_CONFIG%SM%NRFTOTAL
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZGP2DSPP(JL,JRF)=PPERT%PGP2DSPP(JL,1,JRF)
  ENDDO
ENDDO
   
ZFSD(:,:)=0.0_JPRB

! ------------------------------------------------------------------------------
!
!*       2. Call cloud microphysics and convection-cloud interaction
!
! ------------------------------------------------------------------------------
CALL CLOUDSC &
  & (YDECLDP,YDECUMF,YDEPHLI, YDERAD, YDEPHY,  YDVDF, YDSPP_CONFIG,&
  & KDIM%KIDIA,    KDIM%KFDIA,    KDIM%KLON,    KDIM%KLEV, &
  & TSPHY,&
  & STATE%T, STATE%Q, &
  & TENDENCY_CML%T, TENDENCY_CML%Q, TENDENCY_CML%A, TENDENCY_CML%CLD, &
  & TENDENCY_LOC%T, TENDENCY_LOC%Q, TENDENCY_LOC%A, TENDENCY_LOC%CLD, &
  & TENDENCY_VDF%A, TENDENCY_VDF%CLD(:,:,NCLDQL), TENDENCY_VDF%CLD(:,:,NCLDQI),&
  & TENDENCY_DYN%A, TENDENCY_DYN%CLD(:,:,NCLDQL), TENDENCY_DYN%CLD(:,:,NCLDQI),&
  & TENDENCY_DYN%CLD(:,:,NCLDQR), TENDENCY_DYN%CLD(:,:,NCLDQS),&
  & PRAD%PHRSW,    PRAD%PHRLW,&
  & PAUX%PVERVEL,  PAUX%PRSF1,    PAUX%PRS1,&
  & PSURF%PSD_VF(:,YSD_VF%YLSM%MP),PAUX%PGAW,     LLKEYS%LLCUM, &
  & PDIAG%ICTOP,   PDIAG%ITYPE,    PDIAG%IPBLTYPE,PDIAG%ZEIS,&
  & PDIAG%ZLU,     PDIAG%ZLUDE,    PDIAG%ZLUDELI, PDIAG%ZSNDE,   PDIAG%PMFU,     PDIAG%PMFD, ZGP2DSPP,&
  & LDSLPHY, &
  & STATE%A, &
  & STATE%CLD, &
!-- arrays for aerosol-cloud interactions
!!-- aerosol climatology     & ZQAER, 6, &
  & AUXL%ZLCRIT_AER,AUXL%ZICRIT_AER, &     
  & AUXL%ZRE_ICE, &     
  & AUXL%ZCCN,     AUXL%ZNICE,&
!--
  & PDIAG%PCOVPTOT,ZFSD,AUXL%ZRAINFRAC_TOPRFZ,&
  & FLUX%PFCSQL,   FLUX%PFCSQN ,  FLUX%PFCQNNG,  FLUX%PFCQLNG,&
  & FLUX%PFSQRF,   FLUX%PFSQSF ,  FLUX%PFCQRNG,  FLUX%PFCQSNG,&
  & FLUX%PFSQLTUR, FLUX%PFSQITUR , &
  & FLUX%PFPLSL,   FLUX%PFPLSN,   FLUX%PFHPSL,   FLUX%PFHPSN,&
  & PSURF%PSD_XA, KDIM%KFLDX)  

IF(YDEPHY%LRAD_CLOUD_INHOMOG) THEN
  DO JK=1,KDIM%KLEV
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FSD%P(JL,JK) = ZFSD(JL,JK)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUD_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUD_LAYER
