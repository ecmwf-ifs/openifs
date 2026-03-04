! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADFLUX_LAYER(YDSURF,YDMODEL,YDSPP_CONFIG, &
  & KDIM,PAUX,STATE_T0,STATE_TMP,PSURF, SURFL, PRAD,AUXL,FLUX,TENDENCY_LOC, PPERT)

!**** *RADFLUX_LAYER* - Layer routine calling radiative fluxes computation

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ====
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! state_t0  : Derived variable for actual time model state
! state_tmp : Derived variable for updated model state
! PSURF    : Derived variables for general surface quantities

!     ==== Input/output ====
! PRAD     : Derived variables for variables used in radiation scheme
! AUXL     : Derived variables for local auxiliary quantities
! FLUX     : Derived variables for fluxes

!     ==== Output ====
! tendency_loc  :  Derived variables for process tendencies


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
!      Original : 2012-12-04  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      R J Hogan  May-Oct 2014: Pass through surface dn LW, partial derivative,
!                 albedo and PMU0M (cosine of sun angle at last radiation call)
!      M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!      M. Leutbecher (Oct 2020) SPP abstraction
!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM ,   JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_RAD_TYPE, AUX_TYPE, &
   & AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, FLUX_TYPE, PERTURB_TYPE, &
   & SURF_AND_MORE_LOCAL_TYPE
USE SPP_MOD   , ONLY : TSPP_CONFIG
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                    , INTENT(INOUT) :: YDSURF
TYPE(MODEL)                    , INTENT(INOUT) :: YDMODEL
TYPE(TSPP_CONFIG)              , INTENT (IN)   :: YDSPP_CONFIG
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE_T0
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE_TMP
TYPE (SURF_AND_MORE_TYPE)      , INTENT (IN)   :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT (IN)   :: SURFL
TYPE (AUX_RAD_TYPE)            , INTENT(INOUT) :: PRAD
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
TYPE (PERTURB_TYPE)            , INTENT (IN)   :: PPERT
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL, JRF


REAL(KIND=JPRB)    :: ZGP2DSPP(KDIM%KLON, YDSPP_CONFIG%SM%NRFTOTAL)  !SPP: for local 2D random patterns
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "radheatn.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RADFLUX_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, &
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(YSP_RR=>YDSURF%YSP_RR,LERADIMPL=>YDEPHY%LERADIMPL)

!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL RADIATION

! NOTE: LW emissivity PSURF%PEMIS is passed to RADHEATN, and if
! LApproxLwUpdate=TRUE this is then used to update the LW fluxes
! according to the local skin temperature.  However, this was not the
! longwave emissivity actually used by the radiation scheme, which
! constructs emissivity from the contributions from the surface tiles
! and includes separate window and non-window values.  It would be
! desirable to use an emissivity that returns the same upwelling flux
! as the radiation scheme when provided the same skin temperature. To
! attempt this, one could compute this emissivity as a function of the
! surface net longwave (PEMTD), the downwelling longwave (PSRLWD) and
! the skin temperature (PEDRO) at the last call to the radiation
! scheme with:
!    black_body_lw_net = PRAD%PSRLWD(JL) - RSIGMA*PRAD%PEDRO(JL)**4
!    emissivity(JL)    = PRAD%PEMTD(JL,KDIM%KLEV+1)/black_body_lw_net
! However, when this is tried, due to the fact that some of these
! fields have been interpolated back from the radiation grid, the
! resulting emissivity can be a bit noisy (including being outside the
! range 0-1), which is why this is not done at present.

TENDENCY_LOC%T(:,:) = 0._JPRB

!SPP: unpack random pattern fields
DO JRF=1, YDSPP_CONFIG%SM%NRFTOTAL
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZGP2DSPP(JL,JRF)=PPERT%PGP2DSPP(JL,1,JRF)
  ENDDO
ENDDO

IF (LERADIMPL) THEN
CALL RADHEATN &
 & (  YDERAD,YDERDI,YDPHY2,YDPHY3, YDSPP_CONFIG, KDIM%KIDIA  , KDIM%KFDIA  , KDIM%KLON   , KDIM%KLEV, &
 & PAUX%PRS1 ,&
 & PSURF%PEMIS , SURFL%ZSPECTRALEMISS, PRAD%PEMTD, PAUX%PMU0, PAUX%PMU0M,&
 & STATE_TMP%Q,&
 & TENDENCY_LOC%T , PRAD%PTRSW  , AUXL%ZTRSOD , PRAD%PSRLWD, PRAD%PSRSWDC, PRAD%PSRLWDC, &
 & PSURF%PSP_RR(:,YSP_RR%YT%MP9), &
 & PRAD%PHRSW  , PRAD%PHRLW  , PRAD%PHRSC  , PRAD%PHRLC  ,&
 & FLUX%PFRSO  , FLUX%PFRTH  , FLUX%PFRSOD , FLUX%PFRSODC, FLUX%PFRTHD, FLUX%PFRTHDC,&
 & PRAD%PEMTC  , PRAD%PTRSC  , FLUX%PFRSOC , FLUX%PFRTHC , PRAD%PINCSR,&
 & AUXL%ZSUDU  , PRAD%PISUND , PRAD%PDSRP,&
 & PRAD%PSRSWUVB  , PRAD%PSRSWPAR  , PRAD%PSRSWPARC , PRAD%PSRSWTINC , PRAD%PSRFDIR , PRAD%PSRCDIR, &
 & FLUX%PUVDF  , FLUX%PPARF  , FLUX%PPARCF , PRAD%PTINCF , PRAD%PFDIR  , PRAD%PCDIR,  &
 & PRAD%DERIVATIVELW,  &
 & PNEB=PRAD%PNEB   , PAP=PAUX%PRSF1  , PGP2DSPP=ZGP2DSPP)  
ELSE
CALL RADHEATN &
 & (  YDERAD,YDERDI,YDPHY2,YDPHY3, YDSPP_CONFIG, KDIM%KIDIA  , KDIM%KFDIA  , KDIM%KLON   , KDIM%KLEV, &
 & PAUX%PAPRS ,&
 & PSURF%PEMIS , SURFL%ZSPECTRALEMISS, PRAD%PEMTD, PAUX%PMU0, PAUX%PMU0M,&
 & STATE_T0%Q,&
 & TENDENCY_LOC%T , PRAD%PTRSW  , AUXL%ZTRSOD , PRAD%PSRLWD, PRAD%PSRSWDC, PRAD%PSRLWDC, &
 & PSURF%PSP_RR(:,YSP_RR%YT%MP9), &
 & PRAD%PHRSW  , PRAD%PHRLW  , PRAD%PHRSC  , PRAD%PHRLC  ,&
 & FLUX%PFRSO  , FLUX%PFRTH  , FLUX%PFRSOD , FLUX%PFRSODC, FLUX%PFRTHD, FLUX%PFRTHDC,&
 & PRAD%PEMTC  , PRAD%PTRSC  , FLUX%PFRSOC , FLUX%PFRTHC , PRAD%PINCSR,&
 & AUXL%ZSUDU  , PRAD%PISUND , PRAD%PDSRP,&
 & PRAD%PSRSWUVB  , PRAD%PSRSWPAR  , PRAD%PSRSWPARC , PRAD%PSRSWTINC , PRAD%PSRFDIR , PRAD%PSRCDIR, &
 & FLUX%PUVDF  , FLUX%PPARF  , FLUX%PPARCF , PRAD%PTINCF , PRAD%PFDIR  , PRAD%PCDIR,  &
 & PRAD%DERIVATIVELW,  &
 & PNEB=PRAD%PNEB   , PAP=PAUX%PAPRSF  , PGP2DSPP=ZGP2DSPP)  
ENDIF

!       Clear sky fluxes at top and surface
DO JL=KDIM%KIDIA,KDIM%KFDIA
  FLUX%PFRTHC(JL,0)=AUXL%ZCEMTR(JL,0)
  FLUX%PFRTHC(JL,1)=AUXL%ZCEMTR(JL,1)
  FLUX%PFRSOC(JL,0)=PRAD%PINCSR(JL) * AUXL%ZCTRSO(JL,0)
  FLUX%PFRSOC(JL,1)=PRAD%PINCSR(JL) * AUXL%ZCTRSO(JL,1)
  PRAD%PTINCF(JL)  =PRAD%PINCSR(JL)
ENDDO


!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADFLUX_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE RADFLUX_LAYER
