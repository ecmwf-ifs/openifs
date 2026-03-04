! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CHEMINI_LAYER(YDSURF, &
  & YDMODEL,KDIM,PAUX,STATE,PSURF,SURFL,GEMSL)

!**** *CHEMINI_LAYER* - Layer routine calling code to manipulate surface fluxes (injection height, diurnal cycel, deposition)
!                   

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! state    : Derived variable for model state 
! PGFL   - GFL
! PTENGFL - tendency of GFL (GFL time T1)


!     ==== Input/output ====
! PSURF    : Derived variables for general surface quantities
! SURFL    : Derived variables for local surface quantities
! GEMSL    : Derived variable for local GEMS quantities.



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
!      Original : 2014-1-03 ! copy of aerini_layer 

!     MODIFICATIONS.
!     --------------
!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, GEMS_LOCAL_TYPE
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "chem_initflux.intfb.h"
#include "liftemis.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHEMINI_LAYER',0,ZHOOK_HANDLE)


!     ------------------------------------------------------------------
  CALL CHEM_INITFLUX(YDSURF,YDMODEL,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLEV,KDIM%KLON,GEMSL%ITRAC,PSURF%PSD_VF,GEMSL%IAERO, &
 &                   GEMSL%ICHEM,GEMSL%ZCEN,GEMSL%ZTENC,GEMSL%ZCFLX,GEMSL%ZDDVLC, &
 &                   GEMSL%ZSCAV,GEMSL%ZCHEMDV,KDIM%KFLDX,KDIM%KLEVX,&
 &                   PAUX%PGELAM,PAUX%PGELAT,STATE%T(:,KDIM%KLEV), &
 &                   PAUX%PRSF1(:,KDIM%KLEV),PAUX%PDELP,PAUX%PAPHIF,PSURF%PSD_VF(KDIM%KIDIA:KDIM%KFDIA,YDSURF%YSD_VF%YLSM%MP), &
 &                   PSURF%PSD_XA)

 CALL LIFTEMIS(YDMODEL%YRML_GCONF,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,GEMSL%ITRAC,GEMSL%ICHEM,KDIM%KSTGLO, &
 &             PAUX%PDELP,PAUX%PAPHIF,GEMSL%ZCFLX,GEMSL%ZTENC,KDIM%KFLDX,KDIM%KLEVX,PSURF%PSD_XA)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEMINI_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CHEMINI_LAYER
