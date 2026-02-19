! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SURFBC_LAYER(YDSURF, YDEPHY,KDIM, PAUX, &
  & LLKEYS,PSURF,SURFL)

!**** *SURFBC_LAYER* - Layer routine calling SURFBC scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities

!     ==== Input/output ====
! LLKEYS   : Derived variable for local switches
! PSURF    : Derived variable for general surface quantities
! SURFL    : Derived variable for local surface quantities


!        IMPLICIT ARGUMENTS :   NONE
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
!      Original : 2012-11-29  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     A. Agusti-Panareda: 2021-06-21 Pass CO2 photosynthesis type
!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  , ONLY : JPRB,JPIM
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY    , ONLY : TEPHY
USE YOMPHYDER , ONLY : DIMENSION_TYPE, AUX_TYPE, KEYS_LOCAL_TYPE, &
   &                   SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                    , INTENT(INOUT) :: YDSURF
TYPE(TEPHY)                    , INTENT(INOUT) :: YDEPHY
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (KEYS_LOCAL_TYPE)         , INTENT(INOUT) :: LLKEYS
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "surfbc.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFBC_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSURF=>YDEPHY%YSURF, &
 & YSD_VF=>YDSURF%YSD_VF, YSP_SL=>YDSURF%YSP_SL)

!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL SURFBC

CALL SURFBC (YDSURF=YSURF,KIDIA=KDIM%KIDIA,KFDIA=KDIM%KFDIA,KLON=KDIM%KLON,KTILES=KDIM%KTILES,KLEVSN=KDIM%KLEVSN,&
 & PTVL=PSURF%PSD_VF(:,YSD_VF%YTVL%MP),PCO2TYP=PSURF%PSD_VF(:,YSD_VF%YCO2TYP%MP),PTVH=PSURF%PSD_VF(:,YSD_VF%YTVH%MP), &
 & PSOTY=PSURF%PSD_VF(:,YSD_VF%YSOTY%MP),PSDOR=PSURF%PHSTD,PCVLC=PSURF%PSD_VF(:,YSD_VF%YCVL%MP), &
 & PCVHC=PSURF%PSD_VF(:,YSD_VF%YCVH%MP),PCURC=PSURF%PSD_VF(:,YSD_VF%YCUR%MP),&
 & PLAILC=PSURF%PSD_VF(:,YSD_VF%YLAIL%MP),PLAIHC=PSURF%PSD_VF(:,YSD_VF%YLAIH%MP),PLAILI=SURFL%ZLAILI,PLAIHI=SURFL%ZLAIHI,&
 & PLSM=PSURF%PSD_VF(:,YSD_VF%YLSM%MP),PCI=PSURF%PSD_VF(:,YSD_VF%YCI%MP),&
 & PCLAKE=PSURF%PSD_VF(:,YSD_VF%YCLK%MP),PHLICE=PSURF%PSP_SL(:,YSP_SL%YLICD%MP),&
 & PGEMU=PAUX%PGEMU,PSNM1M=SURFL%ZSNM1M,PWLM1M=SURFL%ZWLM1M,PRSNM1M=SURFL%ZRSNM1M,&
 & LDLAND=LLKEYS%LLLAND,LDSICE=LLKEYS%LLSICE,LDLAKE=LLKEYS%LLLAKE,LDNH=LLKEYS%LLNH,LDOCN_KPP=LLKEYS%LLOCN_KPP,&
 & KTVL=PSURF%ITVL,KCO2TYP=PSURF%ICO2TYP,KTVH=PSURF%ITVH,KSOTY=PSURF%ISOTY,&
 & PCVL=PSURF%PCVL,PCVH=PSURF%PCVH,PCUR=PSURF%PCUR,PLAIL=PSURF%PLAIL,PLAIH=PSURF%PLAIH,&
 & PWLMX=SURFL%ZWLMX,PFRTI=SURFL%ZFRTI)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFBC_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE SURFBC_LAYER
