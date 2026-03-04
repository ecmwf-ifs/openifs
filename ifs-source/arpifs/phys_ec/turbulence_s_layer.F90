! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TURBULENCE_S_LAYER(YDSURF, &
 ! Input quantities
  & YDMODEL,KDIM, PAUX, STATE, TENDENCY, AUXL, &
 ! Input/Output quantities
  & GEMSL, PSURF, PCGPP, PCREC, PAG, PRECO, SURFL, FLUX, PDIAG, &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *TURBULENCE_S_LAYER* - Layer routine calling simplified turbulence scheme

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
! tendency : D. V. for model tendencies (entering convection) from processes before 
! AUXL     : Derived variables for local auxiliary quantities

!     ==== Input/output ====
! GEMSL        : Derived variables for local GEMS quantities
! PSURF        : Derived variables for general surface quantities
! PCGPP        : CO2 GPP flux adjustment coefficient
! PCREC        : CO2 REC flux adjustment coefficient
! PAG          : CO2 GPP flux
! PRECO        : CO2 REC flux
! SURFL        : Derived variables for local surface quantities
! FLUX         : Derived variable for fluxes
! PDIAG        : Derived variable for diagnostics quantities

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
!      Original : 2012-11-26  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      2013-08-07 M. Janiskova - adding zfrsoti, zevapsnw to vdfmains
!      E. Dutra/G.Arduini Jan 2018, snow multi-layer arrays, use only 1st level here
!      B. Ingleby  2019-01-17  Move q2 into PSURF%PSD_VD
!      S. Massart   March-2021 : Adding BFAS parameters 

!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER          , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   &                            AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, &
   &                            FLUX_TYPE, GEMS_LOCAL_TYPE,AUX_DIAG_TYPE
USE YOECLDP            , ONLY : NCLDQI,NCLDQL

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT (IN)   :: AUXL
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
REAL(KIND=JPRB)                , INTENT(INOUT) :: PCGPP(KDIM%KLON)
REAL(KIND=JPRB)                , INTENT(INOUT) :: PCREC(KDIM%KLON)
REAL(KIND=JPRB)                , INTENT(INOUT) :: PAG(KDIM%KLON)
REAL(KIND=JPRB)                , INTENT(INOUT) :: PRECO(KDIM%KLON)
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "vdfmains.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TURBULENCE_S_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSP_SG=>YDSURF%YSP_SG, YSP_SB=>YDSURF%YSP_SB, YSD_WS=>YDSURF%YSD_WS, &
 & YSP_RR=>YDSURF%YSP_RR, YSD_VF=>YDSURF%YSD_VF, YSD_VD=>YDSURF%YSD_VD)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL VDFMAINS


CALL VDFMAINS (YDMODEL, KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,KDIM%KLEVS,YDSURF%YSP_SGD%NLEVS,KDIM%KTILES,&
  & GEMSL%ITRAC,TSPHY,PSURF%ITVL,PSURF%ITVH, PSURF%PCVL  , PSURF%PCVH,  PSURF%PCUR,&
  & PSURF%PLAIL, PSURF%PLAIH , PSURF%PFWET , PSURF%PSP_SG(:,:,YSP_SG%YF%MP9), PSURF%PSP_SG(:,:,YSP_SG%YR%MP9), SURFL%ZHSDFOR, &
  & STATE%U    , STATE%V    , STATE%T   , STATE%Q  , GEMSL%ZCEN  , PAUX%PAPRS, PAUX%PAPRSF, PAUX%PGEOM1, PAUX%PGEOMH ,&
  & STATE%CLD(:,:,NCLDQI) , STATE%CLD(:,:,NCLDQL),&
  & PSURF%PSP_RR(:,YSP_RR%YT%MP9)   , PSURF%PSP_SB(:,:,YSP_SB%YT%MP9)  , PSURF%PSP_SB(:,:,YSP_SB%YQ%MP9),&
  & FLUX%PFRSO(:,KDIM%KLEV) , FLUX%PFRTH(:,KDIM%KLEV), PSURF%PEMIS,&
  & SURFL%ZTHKICE, SURFL%ZSNTICE, &
  & PSURF%PSP_SG(:,1,YSP_SG%YT%MP9)  ,PSURF%PSP_SB(:,1,YSP_SB%YTL%MP9),&
  ! & PSURF%PSP_SL(:,YSP_SL%YLICD%MP) , SURFL%ZPTLICE , SURFL%ZPTLWML , & !lake passive
  & PSURF%PSD_VF(:,YSD_VF%YSST%MP)  , PSURF%ISOTY , SURFL%ZFRTI , SURFL%ZALBTI , SURFL%ZWLMX,&
  & PSURF%PSD_WS(:,YSD_WS%YCHAR%MP) ,AUXL%ZTSKRAD, GEMSL%ZCFLX,&
  ! OUTPUT 
  & GEMSL%ZAZ0M , GEMSL%ZAZ0H,&
  & PSURF%PSD_VD(:,YSD_VD%Y10U%MP) , PSURF%PSD_VD(:,YSD_VD%Y10V%MP ) , PSURF%PSD_VD(:,YSD_VD%Y2T%MP) , &
  & PSURF%PSD_VD(:,YSD_VD%Y2D%MP) , PSURF%PSD_VD(:,YSD_VD%Y2SH%MP), &
  & PSURF%PSD_VD(:,YSD_VD%Y10NU%MP), PSURF%PSD_VD(:,YSD_VD%Y10NV%MP),&
  & PDIAG%ZINV,&
  & SURFL%ZFRSOTI, SURFL%ZEVAPSNW, &
  ! INPUT TENDENDCIES
  & TENDENCY%T , TENDENCY%Q , TENDENCY%U , TENDENCY%V , &
  ! OUTPUT TENDENDCIES
  & TENDENCY_LOC%T, TENDENCY_LOC%Q, TENDENCY_LOC%U, TENDENCY_LOC%V, GEMSL%ZTENC, PSURF%PTLE1,&
  !-UPDATED FIELDS FOR TILES 
  & PSURF%PUSTRTI ,PSURF%PVSTRTI ,PSURF%PAHFSTI ,PSURF%PEVAPTI ,PSURF%PTSKTI,&
  ! FLUX OUTPUTS 
  & FLUX%PDIFTS  ,FLUX%PDIFTQ  ,FLUX%PSTRTU  ,FLUX%PSTRTV,&
  & PCGPP, PCREC, PAG, PRECO )

!  Zero the tendencies not used in simplified scheme
TENDENCY_LOC%CLD(:,:,NCLDQL) = 0.0_JPRB
TENDENCY_LOC%CLD(:,:,NCLDQI) = 0.0_JPRB
TENDENCY_LOC%A(:,:)          = 0.0_JPRB

! zeroing fluxes being not computed
FLUX%PSTRDU(KDIM%KIDIA:KDIM%KFDIA,:)=0._JPRB
FLUX%PSTRDV(KDIM%KIDIA:KDIM%KFDIA,:)=0._JPRB

! surface flux of u,v momentum due to gwd 
PDIAG%PUSTRG(KDIM%KIDIA:KDIM%KFDIA)=FLUX%PSTRDU(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV)
PDIAG%PVSTRG(KDIM%KIDIA:KDIM%KFDIA)=FLUX%PSTRDV(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV)


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TURBULENCE_S_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE TURBULENCE_S_LAYER
