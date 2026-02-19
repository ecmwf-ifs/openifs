! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TURBULENCE_LAYER(YDSURF, &
 ! Input quantities
  & YDMODEL, KDIM, PAUX, STATE, TENDENCY, PRAD, PPERT, &
 ! Input/Output quantities
  & AUXL, GEMSL, PSURF, SURFL, FLUX, PDIAG, PERTL, LLKEYS, PDDHS,  &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *TURBULENCE_LAYER* - Layer routine calling turbulence scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! state    : Derived variable for model state
! tendency : D. V. for model tendencies (entering convection) from processes before 
! PRAD     : Derived variable for radiative quantites
! PPERT    : Derived variable for incoming perturbations etc... 

!     ==== Input/output ====
! AUXL     : Derived variable for local auxiliary quantities
! GEMSL    : Derived variable for local GEMS quantities
! PSURF    : Derived variable for general surface quantities
! SURFL    : Derived variable for local surface quantities
! FLUX     : Derived variable for fluxes
! PDIAG    : Derived variable for diagnostics quantities
! PERTL    : Derived variable for local perturbation to physics
! LLKEYS   : Derived variable for local keys
! PDDHS    : Derived variable for surface DDH quantities

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
!      A. Agusti-Panareda 2015-06-10 GPP/REC flux adjustment coefficients
!      M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!      E. Dutra 11-Oct-2016: Add support of snow multi-layer fields from here downwards
!      E.Dutra/G.Arduini Jan 2018: Clean temporary snow multi-layer
!      B. Ingleby  2019-01-17  Move q2 into PSURF%PSD_VD
!      A. Beljaars June 2019: Add optional tile postprocessing
!      P. Bechtold, E. Bazile, I. Sandu 13/01/2020 GFL for ARPEGE-TKE scheme
!      M. Leutbecher Oct 2020 SPP abstraction
!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER    , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   &                      AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, &
   &                      FLUX_TYPE, GEMS_LOCAL_TYPE, AUX_DIAG_TYPE,PERTURB_LOCAL_TYPE, &
   &                      KEYS_LOCAL_TYPE, DDH_SURF_TYPE, AUX_RAD_TYPE, PERTURB_TYPE
USE YOECLDP      , ONLY : NCLDQI, NCLDQL
USE YOMCTESSELDIM, ONLY : IKVTYPES,IKDHVCO2S,IKDHFCO2S,IKDHVVEGS,IKDHFVEGS
USE YOMCST       , ONLY : RCPD

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY
TYPE (AUX_RAD_TYPE)            , INTENT (IN)   :: PRAD
TYPE (PERTURB_TYPE)            , INTENT (IN)   :: PPERT
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (PERTURB_LOCAL_TYPE)      , INTENT(INOUT) :: PERTL
TYPE (KEYS_LOCAL_TYPE)         , INTENT(INOUT) :: LLKEYS
TYPE (DDH_SURF_TYPE)           , INTENT(INOUT) :: PDDHS
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL, JRF, JTILE, JJ, IDIAG
CHARACTER(LEN=1)  :: CLCONF
REAL(KIND=JPRB) :: ZRCPD
REAL(KIND=JPRB) :: ZR0VT(KDIM%KLON), ZEVAPMU(KDIM%KLON)
REAL(KIND=JPRB) :: ZAN(KDIM%KLON)                                    !CTESSEL
REAL(KIND=JPRB) :: ZRD(KDIM%KLON)                                    !CTESSEL
REAL(KIND=JPRB) :: ZRSOIL_STR(KDIM%KLON)                             !CTESSEL
REAL(KIND=JPRB) :: ZDHCO2S(KDIM%KLON,IKVTYPES,IKDHVCO2S+IKDHFCO2S)   !CTESSEL
REAL(KIND=JPRB) :: ZGP2DSPP(KDIM%KLON, YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL) !SPP pattern

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB), ALLOCATABLE :: ZEXDIAG(:,:)

!-----------------------------------------------------------------------

#include "vdfouter.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TURBULENCE_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDMCC=>YDMODEL%YRML_AOC%YRMCC,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
 & YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG)
ASSOCIATE(LBUD23=>YDEPHY%LBUD23, LPPTILES=>YDEPHY%LPPTILES,&
 & LNEMOSICOUP=>YDMCC%LNEMOSICOUP, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSD_WS=>YDSURF%YSD_WS, &
 & YSP_RR=>YDSURF%YSP_RR, YSP_SB=>YDSURF%YSP_SB, YSP_SG=>YDSURF%YSP_SG, &
 & YSP_SL=>YDSURF%YSP_SL, TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

!*         0.     SETUP OF LOCAL ARRAYS

CLCONF='F'  !  Not needed. It just has to be anything else from 'T', so why not 'F' :-)

ZRCPD=1.0_JPRB/RCPD

!   Create space for diagnostic array (for possible pp of canopy resitances)
!   Put LAT/LON in field 1 an 2, just in case that it is needed in SURF for debugging.  
IDIAG=12
ALLOCATE( ZEXDIAG(KDIM%KLON,IDIAG) )
DO JL=KDIM%KIDIA,KDIM%KFDIA
  ZEXDIAG(JL,11)=PAUX%PGELAT(JL) ! Latitude
  ZEXDIAG(JL,12)=PAUX%PGELAM(JL) ! Longitude
ENDDO

!   Stochastic variable with 0. mean and sd 1. (for the moment set to 0; 
!   has to be in the argument list and provided by CALLPAR)

DO JRF=1, YDSPP_CONFIG%SM%NRFTOTAL
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZGP2DSPP(JL,JRF)=PPERT%PGP2DSPP(JL,1,JRF)
  ENDDO
ENDDO

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL VDFMAINS

  CALL VDFOUTER (YDMODEL, CLCONF,&
    & KDIM%KIDIA  , KDIM%KFDIA  , KDIM%KLON   , KDIM%KLEV   , KDIM%KLEVS  , KDIM%KTILES , IKVTYPES, IDIAG, &
    & GEMSL%ITRAC  , GEMSL%ICHEM, GEMSL%IAERO, KDIM%KLEVSN , KDIM%KLEVI  , KDIM%KDHVTLS, KDIM%KDHFTLS, KDIM%KDHVTSS, KDIM%KDHFTSS, &
    & KDIM%KDHVTTS, KDIM%KDHFTTS, KDIM%KDHVTIS, KDIM%KDHFTIS, &
    & IKDHVCO2S,IKDHFCO2S,IKDHVVEGS,IKDHFVEGS, &
    & TSPHY  , PSURF%ITVL   , PSURF%ICO2TYP, PSURF%ITVH  , PSURF%PCVL  , PSURF%PCVH   , PSURF%PCUR, &
    & PSURF%PLAIL  , PSURF%PLAIH  , PSURF%PSD_VF(:,YSD_VF%YFWET%MP) , PSURF%PSP_SG(:,:,YSP_SG%YF%MP9), &
    & PSURF%PSP_SG(:,:,YSP_SG%YR%MP9), PAUX%PMU0 , SURFL%ZHSDFOR, &
    !  & ZR0VT  , &
    & state%u     , state%v     , state%T     , state%q     , state%tke ,&
    & state%cld(:,:,NCLDQL) , state%cld(:,:,NCLDQI) , state%a , GEMSL%ZCEN , &
    & PAUX%PAPRS  , PAUX%PAPRSF , PAUX%PGEOM1 , PAUX%PGEOMH , PAUX%PGELAT , PSURF%PSP_RR(:,YSP_RR%YT%MP9) , &
    & PSURF%PSP_SB(:,:,YSP_SB%YT%MP9) , PSURF%PSP_SB(:,:,YSP_SB%YQ%MP9)   , &
    & FLUX%PFRSO(:,KDIM%KLEV)   , FLUX%PFRTH(:,KDIM%KLEV)   , PSURF%PEMIS  , PRAD%PHRLW  , PRAD%PHRSW  , &
    & PSURF%PSP_SG(:,:,YSP_SG%YT%MP9)  , PSURF%PSP_SB(:,:,YSP_SB%YTL%MP9) , &
    & PSURF%PSP_SL(:,YSP_SL%YLICD%MP9) , PSURF%PSP_SL(:,YSP_SL%YLICT%MP9) , &
    & PSURF%PSP_SL(:,YSP_SL%YLMLT%MP9) , &
    & PSURF%PSD_VF(:,YSD_VF%YSST%MP)   , PSURF%ISOTY  , SURFL%ZFRTI  , SURFL%ZALBTI , SURFL%ZWLMX  , &
    & PSURF%PSD_WS(:,YSD_WS%YCHAR%MP)  , PSURF%PSD_WS(:,YSD_WS%YCHARHQ%MP), &
    & PSURF%PSD_VF(:,YSD_VF%YUCUR%MP)  , PSURF%PSD_VF(:,YSD_VF%YVCUR%MP), &
    & PSURF%PSD_WS(:,YSD_WS%YUSTOKES%MP),PSURF%PSD_WS(:,YSD_WS%YVSTOKES%MP), &
    & AUXL%ZTSKRAD, GEMSL%ZCFLX , GEMSL%ZDDVLC, AUXL%ZSOTEU , AUXL%ZSOTEV , AUXL%ZSOBETA, &
    & SURFL%ZTHKICE, SURFL%ZSNTICE, ZGP2DSPP, &
    & LNEMOSICOUP, &
    ! OUTPUT
    & GEMSL%ZAZ0M  , GEMSL%ZAZ0H  , &
    & PDIAG%PVDIS  , PDIAG%PVDISG, PERTL%ZDISSGW, FLUX%PFTLHEV, FLUX%PFTLHSB, FLUX%PFWSB  , &
    & PSURF%PSD_VD(:,YSD_VD%Y10U%MP)  , PSURF%PSD_VD(:,YSD_VD%Y10V%MP )  , PSURF%PSD_VD(:,YSD_VD%Y2T%MP)  , &
    & PSURF%PSD_VD(:,YSD_VD%Y2D%MP)  , PSURF%PSD_VD(:,YSD_VD%Y2SH%MP) , PDIAG%ZINV  , &
    & PSURF%PSD_VD(:,YSD_VD%YBLH%MP) , PDIAG%ZEIS, &
    & PSURF%PSD_VD(:,YSD_VD%Y10NU%MP) , PSURF%PSD_VD(:,YSD_VD%Y10NV%MP) , PSURF%PSD_VD(:,YSD_VD%YZUST%MP) , &
    & SURFL%ZFRSOTI, SURFL%ZEVAPSNW , ZEVAPMU , ZEXDIAG, PDIAG%PI10FG  , PSURF%PSD_VD(:,YSD_VD%Y10FGCV%MP), PDIAG%IPBLTYPE, &
    ! Input tendencies
    & tendency%T, tendency%q, tendency%cld(:,:,NCLDQL), tendency%cld(:,:,NCLDQI), &
    & tendency%a, tendency%u, tendency%v, &
    ! OUTPUT TENDENDCIES
    & tendency_loc%T, tendency_loc%q, tendency_loc%cld(:,:,NCLDQL), tendency_loc%cld(:,:,NCLDQI), &
    & tendency_loc%a, tendency_loc%u, tendency_loc%v, tendency_loc%tke,&
    & GEMSL%ZTENC  , PSURF%PTLE1  , &
    ! UPDATED FIELDS FOR TILES
    & PSURF%PUSTRTI, PSURF%PVSTRTI, PSURF%PAHFSTI, PSURF%PEVAPTI, PSURF%PTSKTI , &
    & SURFL%ZAHFTRTI, &
    !-UPDATED FIELDS FOR VEGETATION TYPES
    & SURFL%ZANDAYVT,SURFL%ZANFMVT,&
    ! FLUX OUTPUTS
    & FLUX%PDIFTS , FLUX%PDIFTQ , FLUX%PDIFTL , FLUX%PDIFTI , FLUX%PSTRTU , FLUX%PSTRTV , FLUX%PTOFDU , FLUX%PTOFDV , &
    & FLUX%PSTRDU , FLUX%PSTRDV , PDIAG%PKH_VDF, PDIAG%PKM_VDF, PDIAG%ZRI, &
    & ZAN    , PSURF%PSD_VD(:,YSD_VD%YIGPP%MP)    , ZRD    , ZRSOIL_STR    , &
    & PSURF%PSD_VD(:,YSD_VD%YIREC%MP), PSURF%PSD_VD(:,YSD_VD%YINEE%MP), PSURF%PSD_VD(:,YSD_VD%YICH4%MP), GEMSL%ZCFLXO, & 
    ! DDH OUTPUTS
    & PDDHS%PDHTLS , PDDHS%PDHTSS , PDDHS%PDHTTS , PDDHS%PDHTIS, &
    & ZDHCO2S, SURFL%ZDHVEGS, &
    ! BFAS ADJUSTMENT COEFFICIENTS FOR CTESSEL CO2 VEGETATION FLUXES
    & PCGPP=PSURF%PSD_VF(:,YSD_VF%YCGPP%MP),PCREC=PSURF%PSD_VF(:,YSD_VF%YCREC%MP))

! key for shallow convection (always on)
DO JL=KDIM%KIDIA,KDIM%KFDIA
  LLKEYS%LLSHCV(JL) = .TRUE.
ENDDO

!   Diagnostics using local arrays
IF(LBUD23) THEN
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,12)=PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,12) &
   &                     +PERTL%ZDISSGW(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)*ZRCPD
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_XA(JL,7,25) = PDIAG%ZINV(JL)
    PSURF%PSD_XA(JL,8,25) = PDIAG%ZEIS(JL)
    PSURF%PSD_XA(JL,9,25) = PDIAG%IPBLTYPE(JL)
    SELECT CASE (PDIAG%IPBLTYPE(JL))
    CASE (0)
    PSURF%PSD_XA(JL,10,25) = PSURF%PSD_XA(JL,10,25) + 1.0
    CASE (1)
    PSURF%PSD_XA(JL,11,25) = PSURF%PSD_XA(JL,11,25) + 1.0
    CASE (2)
    PSURF%PSD_XA(JL,12,25) = PSURF%PSD_XA(JL,12,25) + 1.0
    CASE (3)
    PSURF%PSD_XA(JL,13,25) = PSURF%PSD_XA(JL,13,25) + 1.0
    END SELECT
  ENDDO
ENDIF

! surface flux of u,v momentum due to gwd 
DO JL=KDIM%KIDIA,KDIM%KFDIA
  PDIAG%PUSTRG(JL)=FLUX%PSTRDU(JL,KDIM%KLEV)
  PDIAG%PVSTRG(JL)=FLUX%PSTRDV(JL,KDIM%KLEV)
ENDDO

! Put ZEVAPMU (potential evaporation) into new field for pp
DO JL=KDIM%KIDIA,KDIM%KFDIA
  IF (PSURF%PSD_VF(JL,YSD_VF%YLSM%MP) > 0.5_JPRB) THEN
    PSURF%PEVAPMU(JL) = ZEVAPMU(JL)
  ELSE
    PSURF%PEVAPMU(JL) = 0.0_JPRB
  ENDIF
ENDDO

! Additonal diagnostics + storage of tiled fluxes (KLEV has to be >=137)
IF (LPPTILES) THEN
  IF (KDIM%KLEV < 137) THEN
    CALL ABOR1("TURBULENCE_LAYER: PPTILES requires KLEV >= 137")
  ENDIF

  PSURF%PSD_XA(:,:,1)=0.0_JPRB 

  ! Extra surface diagnostics 
  DO JJ=1,10
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PSURF%PSD_XA(JL,JJ+100,1) = ZEXDIAG(JL,JJ)
    ENDDO
  ENDDO

  ! Lowest model level prognostic variables 
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_XA(JL,130,1) = state%u(JL,KDIM%KLEV)
    PSURF%PSD_XA(JL,131,1) = state%v(JL,KDIM%KLEV)
    PSURF%PSD_XA(JL,132,1) = state%T(JL,KDIM%KLEV)
    PSURF%PSD_XA(JL,133,1) = state%q(JL,KDIM%KLEV)
  ENDDO

  ! Tiled variables 
  DO JTILE=1,KDIM%KTILES
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PSURF%PSD_XA(JL,   JTILE,1) = PSURF%PEVAPTI(JL,JTILE)
      PSURF%PSD_XA(JL,10+JTILE,1) = PSURF%PAHFSTI(JL,JTILE)
      PSURF%PSD_XA(JL,20+JTILE,1) = PSURF%PTSKTI(JL,JTILE)
      PSURF%PSD_XA(JL,30+JTILE,1) = SURFL%ZAHFTRTI(JL,JTILE)
      PSURF%PSD_XA(JL,40+JTILE,1) = SURFL%ZFRSOTI(JL,JTILE)
      PSURF%PSD_XA(JL,50+JTILE,1) = SURFL%ZFRTI(JL,JTILE)
    ENDDO
  ENDDO
ENDIF

! Release memory
DEALLOCATE (ZEXDIAG)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TURBULENCE_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE TURBULENCE_LAYER
