! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLOUD_S_LAYER(&
 ! Input quantities
 &  YDMODEL,KDIM, LDRAIN1D, PAUX, STATE, TENDENCY, &
 ! Input/Output quantities
 &  AUXL, FLUX, PDIAG,&
 ! Output tendencies
 &  TENDENCY_LOC)

!**** *CLOUD_S_LAYER* - Layer routine calling simplified cloud scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDRAIN1D : key to activate 1DRAIN
! state    : Derived variable for model state
! tendency : D. V. for model tendencies (entering cloud) from processes before 
! PAUX     : Derived variables for general auxiliary quantities

!     ==== Input/output ====
! AUXL         : Derived variable for local quantites
! FLUX         : Derived variable for fluxes
! PDIAG        : Derived variable for diagnostics

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
!      F. Vana  8-Apr-2016  small fix
!      F. Vana  14-Sep-2020 arguments update & simplified prognostic scheme 

!-----------------------------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &                    FLUX_TYPE, MASK_GFL_TYPE, AUX_DIAG_LOCAL_TYPE, AUX_DIAG_TYPE
USE YOECLDP   , ONLY : NCLDQR,NCLDQS,NCLDQI,NCLDQL

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)                     ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
LOGICAL                        , INTENT (IN)   :: LDRAIN1D
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFLAG
REAL(KIND=JPRB)    :: ZCONS
REAL(KIND=JPRB)    :: Z_T1(KDIM%KLON,KDIM%KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "satur.intfb.h"
#include "cloudst.intfb.h"
#include "cloudsc2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUD_S_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YDPHNC=>YDMODEL%YRML_PHY_SLIN%YRPHNC)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY,LEPCLD2=>YDPHNC%LEPCLD2, &
 &   LPHYLIN=>YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN)
!     ------------------------------------------------------------------

!*         0.     RECOMPUTE SATURATION TO APPROPRIATE STATE
IFLAG=2
Z_T1(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)=STATE%T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
 & + TSPHY*TENDENCY%T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
CALL SATUR (KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KTDIA, KDIM%KLEV, LPHYLIN, &
 & PAUX%PRSF1, Z_T1, PDIAG%PQSAT, IFLAG) 


!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL SIMPLIFIED CLOUD

IF (LEPCLD2) THEN

  !  PROGNOSTIC SCHEME
  CALL CLOUDSC2 (&
    & YDMODEL%YRML_PHY_SLIN%YREPHLI,YDPHNC, &
    & YDMODEL%YRML_PHY_EC%YRECLD,YDMODEL%YRML_PHY_EC%YRECLDP, &
    & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON , KDIM%KTDIA , KDIM%KLEV, LDRAIN1D, &
    & TSPHY  , &
    & PAUX%PRS1  , PAUX%PRSF1, STATE%Q   , PDIAG%PQSAT , STATE%T, &
    & STATE%CLD(:,:,NCLDQL),STATE%CLD(:,:,NCLDQI), &
    & PDIAG%ZLUDE , PDIAG%ZLU  , PDIAG%PMFU, PDIAG%PMFD, &
    & TENDENCY_LOC%T, TENDENCY%T, TENDENCY_LOC%Q, TENDENCY%Q, &
    & TENDENCY_LOC%CLD(:,:,NCLDQL), TENDENCY%CLD(:,:,NCLDQL), & 
    & TENDENCY_LOC%CLD(:,:,NCLDQI), TENDENCY%CLD(:,:,NCLDQI), &
    & AUXL%ZPATMP ,FLUX%PFPLSL  , FLUX%PFPLSN, &
    & FLUX%PFHPSL, FLUX%PFHPSN, PDIAG%PCOVPTOT )
ELSE

  !  DIAGNOSTIC SCHEME

  CALL CLOUDST &
    & ( YDMODEL%YRML_PHY_SLIN%YREPHLI,YDPHNC, &
    & YDMODEL%YRML_PHY_EC%YRECLD,YDMODEL%YRML_PHY_EC%YRECLDP, &
    & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON , KDIM%KTDIA , KDIM%KLEV, LDRAIN1D, &
    & TSPHY  , &
    & PAUX%PRS1  , PAUX%PRSF1, STATE%Q   , PDIAG%PQSAT , STATE%T, &
    & PDIAG%ZLUDE , PDIAG%ZLU   , &
    & TENDENCY_LOC%T, TENDENCY%T, TENDENCY_LOC%Q, TENDENCY%Q, &
    & AUXL%ZPATMP , AUXL%ZPLSMP, AUXL%ZPLTMP,FLUX%PFPLSL  , FLUX%PFPLSN, &
    & FLUX%PFHPSL, FLUX%PFHPSN, PDIAG%PCOVPTOT )

ENDIF

!*         2.     SECURITY AND WRAP UP

! Turn semi-prog quantities to tendency and zero precipitations
ZCONS=1._JPRB/TSPHY
IF (TENDENCY_LOC%KEYMASK%A ) TENDENCY_LOC%A(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
 &  = ZCONS*(AUXL%ZPATMP(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)                    &
 &             - STATE%A(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV))
IF (.NOT.LEPCLD2 .AND. TENDENCY_LOC%KEYMASK%QL)                         &
 & TENDENCY_LOC%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQL)           &
 &  = ZCONS*(AUXL%ZPLTMP(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)             &
 &            - STATE%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQL))
IF (.NOT.LEPCLD2 .AND. TENDENCY_LOC%KEYMASK%QI)                         &
 & TENDENCY_LOC%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQI)           &
 &  = ZCONS*(AUXL%ZPLSMP(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)             &
 &            - STATE%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQI))
IF (TENDENCY_LOC%KEYMASK%QR) TENDENCY_LOC%CLD(:,:,NCLDQR)=0.0_JPRB
IF (TENDENCY_LOC%KEYMASK%QS) TENDENCY_LOC%CLD(:,:,NCLDQS)=0.0_JPRB

! Some security for diagnostics
PDIAG%PPRECTYPE(KDIM%KIDIA:KDIM%KFDIA) = 0._JPRB
PDIAG%PFZRA    (KDIM%KIDIA:KDIM%KFDIA) = 0._JPRB
!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUD_S_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUD_S_LAYER
