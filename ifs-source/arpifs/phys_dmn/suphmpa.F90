! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUPHMPA(YDGEOMETRY,YDLDDH,YDML_GCONF,YDDYNA,YDML_PHY_MF,KULOUT)

!**** *SUPHMPA*   - Initialize common meso_NH MODD_ used in physics for AROME

!     Purpose.
!     --------
!           Initialize MODD_PARAMETERS, MODD_CST, MODD_CONF,
!           MODD_RAIN_ICE_DESCR, MODD_RAIN_ICE_PARAM, MODD_BUDGET  
!           parameters used in meso_NH Physics and aladin/meso_NH physics
!           interface 

!**   Interface.
!     ----------
!        *CALL* *SUPHMPA(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!       R. Zaaboul
!        Original : 28-Feb-2006

!     Modifications.
!     --------------
!           E. BAZILE : 01-09-2006 : Modified for LCVPPKF.
!        Y. Seity : 28-March-2007 Move chemistry setup under suphmse       
!        O.Riviere: 01/10/2008 removal of call to now obsolete aro_iniapft
!        S.Riette: 24 Aug 2011 add call to AROINI_NEB
!        Y.Seity: 9 Feb 2014 : add autoconversion setup (*CRIAUT*)
!        Y.Seity: 12 Nov  2014 : add test on NGFL_EZDIAG
!        S. Riette (Jan 2015): new ICE3 and ICE4 parameters with new aroini_micro interface
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMDYNA  , ONLY  : TDYNA
USE YOMLDDH   , ONLY : TLDDH
USE YOMCT0 ,ONLY : LELAM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TLDDH)       ,INTENT(INOUT) :: YDLDDH
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TDYNA), INTENT(IN)          :: YDDYNA
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM)   :: ILON
INTEGER(KIND=JPIM)   :: ILEV
INTEGER(KIND=JPIM)   :: ISV
INTEGER(KIND=JPIM)   :: IRR

REAL(KIND=JPRB) :: ZTSTEP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-- Next five strings are kept in modules in MNH, not in AROME !!!
!   ( here they are hardcoded for use in AROINI_BUDGET )

CHARACTER (LEN=4)   :: CLUOUT     ! name of output listing
CHARACTER (LEN=4)   :: CLTURB     ! Kind of turbulence parameterization
CHARACTER (LEN=4)   :: CLTURBDIM  ! dimensionality of the turbulence scheme
CHARACTER (LEN=4)   :: CLRAD      ! Kind of radiation scheme
CHARACTER (LEN=4)   :: CLDCONV    ! Kind of deep convection scheme
CHARACTER (LEN=4)   :: CLSCONV    ! Kind of shallow convection scheme
CHARACTER (LEN=4)   :: CLCLOUD    ! Kind of microphysical scheme
CHARACTER (LEN=4)   :: CLMET_ADV_SCHEME ! type of advection scheme for meteorological scalar variables
CHARACTER (LEN=4)   :: CLSV_ADV_SCHEME  ! type of advection scheme for tracer scalar variables

LOGICAL :: LLNOTMAP
!     ------------------------------------------------------------------

#include "sucvmnh.intfb.h"
#include "aroini_cstmnh.h"
#include "aroini_micro.h"
#include "aro_subudget.h"
#include "aroini_budget.h"
#include "aroini_turb.h"
#include "abor1.intfb.h"
#include "aroini_mfshal.h"

#include "aroini_micro_lima.h"

IF (LHOOK) CALL DR_HOOK('SUPHMPA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDPHY=>YDML_PHY_MF%YRPHY,YDRIP=>YDML_GCONF%YRRIP,YDARPHY=>YDML_PHY_MF%YRARPHY,YDPARAR=>YDML_PHY_MF%YRPARAR)
ASSOCIATE(XDETR_LUP=>YDPARAR%XDETR_LUP, XCMF=>YDPARAR%XCMF, &
 & XBDETR=>YDPARAR%XBDETR, XLINI=>YDPARAR%XLINI, XABUO=>YDPARAR%XABUO, &
 & XLAMBDA=>YDPARAR%XLAMBDA, RT0CRIAUTI=>YDPARAR%RT0CRIAUTI, &
 & XKCF_MF=>YDPARAR%XKCF_MF, XALP_PERT=>YDPARAR%XALP_PERT, &
 & NSPLITR=>YDPARAR%NSPLITR, NSPLITG=>YDPARAR%NSPLITG, XSIGMA_MF=>YDPARAR%XSIGMA_MF, XA1=>YDPARAR%XA1, &
 & CMICRO=>YDPARAR%CMICRO, XENTR_DRY=>YDPARAR%XENTR_DRY, &
 & XENTR_MF=>YDPARAR%XENTR_MF, CSEDIM=>YDPARAR%CSEDIM, NSV=>YDPARAR%NSV, &
 & XFRAC_UP_MAX=>YDPARAR%XFRAC_UP_MAX, XB=>YDPARAR%XB, XC=>YDPARAR%XC, &
 & XTAUSIGMF=>YDPARAR%XTAUSIGMF, RCRIAUTC=>YDPARAR%RCRIAUTC, &
 & XDETR_DRY=>YDPARAR%XDETR_DRY, XR=>YDPARAR%XR, RCRIAUTI=>YDPARAR%RCRIAUTI, &
 & LOWARM=>YDPARAR%LOWARM, XBENTR=>YDPARAR%XBENTR, XBETA1=>YDPARAR%XBETA1, &
 & LAROBU_ENABLE=>YDPARAR%LAROBU_ENABLE, LCRIAUTI=>YDPARAR%LCRIAUTI, &
 & XKRC_MF=>YDPARAR%XKRC_MF, XALPHA_MF=>YDPARAR%XALPHA_MF, &
 & XPRES_UV=>YDPARAR%XPRES_UV, NRR=>YDPARAR%NRR, XCRAD_MF=>YDPARAR%XCRAD_MF, &
 & CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, LHARATU=>YDPARAR%LHARATU, &
 & LMPA=>YDARPHY%LMPA, LKFBCONV=>YDARPHY%LKFBCONV, LMFSHAL=>YDARPHY%LMFSHAL, &
 & LMICRO=>YDARPHY%LMICRO, LTURB=>YDARPHY%LTURB, &
 & NPROMA=>YDDIM%NPROMA, &
 & LEDKF=>YDPHY%LEDKF, LCVPPKF=>YDPHY%LCVPPKF, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LSDDH=>YDLDDH%LSDDH, TSTEP=>YDRIP%TSTEP, &
 & XTSTEP_TS=>YDPARAR%XTSTEP_TS, CSNOWRIMING=>YDPARAR%CSNOWRIMING, XMRSTEP=>YDPARAR%XMRSTEP,&
 & NMAXITER_MICRO=>YDPARAR%NMAXITER_MICRO, LFEEDBACKT=>YDPARAR%LFEEDBACKT, LEVLIMIT=>YDPARAR%LEVLIMIT, &
 & LNULLWETG=>YDPARAR%LNULLWETG, LWETGPOST=>YDPARAR%LWETGPOST, LNULLWETH=>YDPARAR%LNULLWETH,&
 & LWETHPOST=>YDPARAR%LWETHPOST, XFRACM90=>YDPARAR%XFRACM90, LCONVHG=>YDPARAR%LCONVHG, &
 & CSUBG_RC_RR_ACCR=>YDPARAR%CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP=>YDPARAR%CSUBG_RR_EVAP, &
 & CSUBG_PR_PDF=>YDPARAR%CSUBG_PR_PDF, LCRFLIMIT=>YDPARAR%LCRFLIMIT, &
 & CFRAC_ICE_ADJUST=>YDPARAR%CFRAC_ICE_ADJUST, XSPLIT_MAXCFL=>YDPARAR%XSPLIT_MAXCFL,&
 & CFRAC_ICE_SHALLOW_MF=>YDPARAR%CFRAC_ICE_SHALLOW_MF, LSEDIM_AFTER=>YDPARAR%LSEDIM_AFTER)
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!       1. Initialisation of MesoNH constantes

IF (LELAM) THEN
  LLNOTMAP=.NOT.YDGEOMETRY%YREGEO%LMAP
ELSE
  LLNOTMAP=.TRUE.
ENDIF
CALL AROINI_CSTMNH (KULOUT,YDDYNA%LTWOTL,LLNOTMAP)

!       2. Initialisation for microphysics scheme

ZTSTEP=TSTEP
IF (CMICRO == 'LIMA') THEN
   CALL AROINI_MICRO (KULOUT,ZTSTEP,LOWARM,'ICE3',NSPLITR,CSEDIM,LCRIAUTI,&
        &             RCRIAUTI,RT0CRIAUTI,RCRIAUTC, XTSTEP_TS, CSNOWRIMING, XMRSTEP,&
 &                    NMAXITER_MICRO, LFEEDBACKT, LEVLIMIT, LNULLWETG, LWETGPOST, LNULLWETH,&
 &                    LWETHPOST, XFRACM90, LCONVHG, CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP, CSUBG_PR_PDF,&
 &                    LCRFLIMIT, CFRAC_ICE_ADJUST, XSPLIT_MAXCFL,&
 &                    CFRAC_ICE_SHALLOW_MF, LSEDIM_AFTER)
   CALL AROINI_MICRO_LIMA (KULOUT,4,ZTSTEP,LOWARM,CMICRO,NSPLITR,NSPLITG,CSEDIM,LCRIAUTI,&
        &                 RCRIAUTI,RT0CRIAUTI,RCRIAUTC)
ELSE
CALL AROINI_MICRO (KULOUT,ZTSTEP,LOWARM,CMICRO,NSPLITR,CSEDIM,LCRIAUTI,&
 &                 RCRIAUTI,RT0CRIAUTI,RCRIAUTC,XTSTEP_TS, CSNOWRIMING, XMRSTEP,&
 &                 NMAXITER_MICRO, LFEEDBACKT, LEVLIMIT, LNULLWETG, LWETGPOST, LNULLWETH,&
 &                 LWETHPOST, XFRACM90, LCONVHG, CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP, CSUBG_PR_PDF,&
 &                 LCRFLIMIT, CFRAC_ICE_ADJUST, XSPLIT_MAXCFL,&
 &                 CFRAC_ICE_SHALLOW_MF, LSEDIM_AFTER)
 ENDIF

!       3. Initialisation of Budget

ILON = NPROMA
ILEV = NFLEVG
ISV  = NSV
IRR  = NRR

IF (LMPA.AND.LSDDH) THEN
LAROBU_ENABLE=.TRUE.
CALL ARO_SUBUDGET(ILON,ILEV,ZTSTEP)
ENDIF
 
!! WARNING:
!!    Values of CLRAD,CLDCONV,CLTURB,CLTURBDIM, CLCLOUD
!!    are hard coded here. In MNH they are saved
!!    in modules.
 
CLUOUT     = 'AIBUO'
CLRAD      = 'NONE'
CLSCONV    = CMF_UPDRAFT
CLMET_ADV_SCHEME='XXXX'
CLSV_ADV_SCHEME='XXXX'

IF(LKFBCONV.OR.LCVPPKF) THEN
  CLDCONV    = 'KAFR'
ELSE
  CLDCONV    = 'NONE'
ENDIF
 
IF(LTURB) THEN
  CLTURB    = 'XXXX'
  CLTURBDIM = '1DIM'
ELSE
  CLTURB    = 'NONE'
  CLTURBDIM = 'NONE'
ENDIF
 
IF(LMICRO) THEN
  CLCLOUD = CMICRO
ELSE
  CLCLOUD = 'NONE'
ENDIF
 
IF (LMPA.AND.LSDDH) THEN
CALL AROINI_BUDGET(LAROBU_ENABLE,CLUOUT,KULOUT,ZTSTEP,ISV,IRR,&
                   & CLRAD,CLDCONV,CLSCONV,CLTURB,CLTURBDIM,&
                   & CLCLOUD,CLMET_ADV_SCHEME, CLSV_ADV_SCHEME)


ENDIF


!       4. Initialisation of Turbulence scheme

CALL AROINI_TURB(XLINI,LHARATU)

!       5. Initialisation of Mass Flux Shallow convection scheme

IF(LMFSHAL.OR.LEDKF) CALL AROINI_MFSHAL(XALP_PERT,XABUO,XBENTR,XBDETR,XCMF,XENTR_MF,XCRAD_MF,XENTR_DRY,&
 &          XDETR_DRY,XDETR_LUP,XKCF_MF,XKRC_MF,XTAUSIGMF,XPRES_UV,XFRAC_UP_MAX,&
 &          XALPHA_MF,XSIGMA_MF,XA1,XB,XC,XBETA1,XR,XLAMBDA)

IF (LMFSHAL.AND.YDML_GCONF%YGFL%NGFL_EZDIAG < 3) THEN
  CALL ABOR1 ("With LMFSHAL NGFL_EZDIAG should be >= 3 !")
ENDIF

!       6. Initialisation of Convection scheme

IF(LKFBCONV.OR.LCVPPKF) THEN
  CALL SUCVMNH(YDML_PHY_MF,KULOUT)
ENDIF

!       7. Initialisation of nebulosity computation

CALL AROINI_NEB

! -----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHMPA',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHMPA
