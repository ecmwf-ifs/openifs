! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE SUINIF(YDGEOMETRY,YDGFL,YDSURF,YDSP,YDCFU,YDXFU,YDMODEL,KFILE,LDRDGRIDSP,YDMCUF,KINITMONTH)

!**** *SUINIF*  - Routine to initialize the fields of the model.

!     Purpose.
!     --------
!           Initialize the fields of the model.

!**   Interface.
!     ----------
!        *CALL* *SUINIF(...)*

!        Explicit arguments :
!        --------------------
!        KFILE : an indicator for which  file is to be read
!               KFILE = 0 the CFNISH file is read
!               KFILE = 1 the CFNGSH file is read
!               KFILE = 2 the CFNRF  file is read
!               KFILE = 3 the CFNBS file is read (bias, for recursive dfi)
!               KFILE = 4 the CFNGG  file is read
!               KFILE = 6 the CFNFGI file is read
!               KFILE = 7 the CFNANI file is read
!               KFILE = 8 the CFNPANSH file is read
!               KFILE =12 the CFNBGHRSH/GG file is read
!        LDRDGRIDSP: T (resp. F): GMV, GMVS and spectral GFL read as grid-point
!               fields (resp. spectral fields) on files.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!        Mats Hamrud   *ECMWF*
!        Original : 91-11-14

!     Modifications.
!     --------------
!        N.Wedi/K.Yessad (Jan 2008): different dev for NH model and PC scheme
!        E.Holm        13-Nov-2008 Add high resolution background=12, used
!                                  by rdfpinc for adding increments
!        N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!        R. El Khatib : 09-99-2008 : fix use of LNHDYN+LSPRT
!        R. El Khatib  26-Nov-2008  Cleaning
!        Y.Takaya      29-Jan-2009 Ocean mixed layer fields
!        K. Yessad (July 2014): Move some variables.
!        R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!        O. Marsden      Aug 2016 Call cleanups due to removal of SPA3
!        K. Yessad (Feb 2018): remove deep-layer formulations.
!     -----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGFL             , ONLY : TGFL
USE YOMCFU             , ONLY : TCFU
USE YOMCST             , ONLY : RTT
USE YOMXFU             , ONLY : TXFU
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : NCONF, LFBDAP, LELAM, LIFSMIN, LECMWF
USE YOMLUN             , ONLY : NULOUT
USE TRAJECTORY_MOD     , ONLY : LTRAJHR
!USE SPECTRAL_FIELDS_MOD, ONLY: SPECTRAL_FIELD, ASSIGNMENT(=)
USE SPECTRAL_FIELDS_MOD
USE YOMMCUF            , ONLY : TMCUF
USE MPL_MODULE         , ONLY : MPL_BARRIER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for call to SUSPEC
TYPE(TGFL)          , INTENT(INOUT) :: YDGFL
TYPE(TSURF)         , INTENT(INOUT) :: YDSURF
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP
TYPE(TCFU)          , INTENT(INOUT) :: YDCFU
TYPE(TXFU)          , INTENT(INOUT) :: YDXFU
TYPE(MODEL)         , INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM)  , INTENT(IN)    :: KFILE
LOGICAL, OPTIONAL   , INTENT(IN)    :: LDRDGRIDSP
TYPE(TMCUF)         , INTENT(INOUT), OPTIONAL :: YDMCUF
INTEGER(KIND=JPIM), INTENT(IN),    OPTIONAL :: KINITMONTH

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFILE
INTEGER(KIND=JPIM) :: J, JSTGLO, ICEND, IBL, JLEV, ILEVOML 
LOGICAL :: LLINOR, LLSPOR
LOGICAL :: LLRDGRIDSP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "especrt.intfb.h"
#include "specrt.intfb.h"
#include "spnh_conv_nhvar.intfb.h"
#include "sugrcfu.intfb.h"
#include "sugridf.intfb.h"
#include "sugridu.intfb.h"
#include "sugrido.intfb.h"
#include "sugrxfu.intfb.h"
#include "suspec.intfb.h"
#include "spnorm.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUINIF',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDPHNC=>YDMODEL%YRML_PHY_SLIN%YRPHNC)
ASSOCIATE(NUMGPFLDS=>YGFL%NUMGPFLDS, &
 & LCUMFU=>YDCFU%LCUMFU, LREACFU=>YDCFU%LREACFU, &
 & LEOCML=>YDEPHY%LEOCML, &
 & LOCMLTKE=>YDEPHY%LOCMLTKE, &
 & NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & SD_VF=>YDSURF%SD_VF, YSD_VF=>YDSURF%YSD_VF, &
 & YSP_OMD=>YDSURF%YSP_OMD, &
 & SP_OM=>YDSURF%SP_OM, YSP_OM=>YDSURF%YSP_OM, &
 & SD_V3=>YDSURF%SD_V3, YSD_V3=>YDSURF%YSD_V3, &
 & LENCLD2=>YDPHNC%LENCLD2, LEPCLD2=>YDPHNC%LEPCLD2, &
 & LREAXFU=>YDXFU%LREAXFU, LXFU=>YDXFU%LXFU, &
 & LREASUR=>YDPHY%LREASUR)
!     ------------------------------------------------------------------

!*       0.    PRELIMINARY SPECTRAL FIT OF PRE-PROCESSED DATA
!              ----------------------------------------------

CALL GSTATS(14,0)

IF (PRESENT(LDRDGRIDSP)) THEN
  LLRDGRIDSP=LDRDGRIDSP
ELSE
  LLRDGRIDSP=.FALSE.
ENDIF

!*       1.    INITIALIZE SPECTRAL FIELDS
!              ---------------------------

LLINOR=(NCONF /= 131. .OR. .NOT.LECMWF).AND.(KFILE /= 6 .AND. KFILE/10 /= 6 .AND. KFILE /= 7)

CALL SUSPEC(YDGEOMETRY,YDGFL,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYN,YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_LBC, &
 &          YDSP,KFILE=KFILE,LDRDGRIDSP=LLRDGRIDSP,LDINOR=LLINOR, LDSPOR=LLSPOR, YDMCUF=YDMCUF)
WRITE(NULOUT,'(A,I3,A)') ' SUINIF: SUSPEC(KFILE=',KFILE,') called '
CALL FLUSH(NULOUT)

CALL MPL_BARRIER(CDSTRING='SUINIF1')

!     ------------------------------------------------------------------

!*       2.    INITIALIZE GRID POINT FIELDS
!              ----------------------------

!*       2.1   SURFACE
IF(LREASUR) THEN
  CALL SUGRIDF(YDGEOMETRY,YDSURF,YDMODEL,KFILE,YDSP%OROG,LLSPOR,KINITMONTH=KINITMONTH)
ENDIF

!*       2.2   OCEAN MIXED LAYER (KPP)
IF(LEOCML) THEN
  CALL SUGRIDO(YDGEOMETRY,YDSURF,YDEPHY,KFILE)
ELSEIF (LOCMLTKE) THEN
  !! the Ocean TKE scheme has so far been implemented without initial conditions
  !! Initialise here with temperature profile given by SST, and no motion (Salinity is not used but initialised to 0.)
  !! and ocean TKE = 0. (see surf/module/oc_mlm_mod.F90 for the actul initialisation of TKE)
  ILEVOML=YSP_OMD%NLEVS

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSTGLO,ICEND,JLEV,J,IBL)
  DO JSTGLO=1,NGPTOT,NPROMA
    ICEND = MIN(NPROMA,NGPTOT-JSTGLO+1)
    IBL   = (JSTGLO-1)/NPROMA+1
    DO JLEV = 1, ILEVOML
      DO J = 1, ICEND
        SP_OM(J,JLEV,YSP_OM%YTO%MP,IBL) = SD_VF(J,YSD_VF%YSST%MP,IBL) - RTT  !!!  YTO is in degree C, YSST in Kelvin !!!
        SP_OM(J,JLEV,YSP_OM%YSO%MP,IBL) = 0.0_JPRB
        SP_OM(J,JLEV,YSP_OM%YUO%MP,IBL) = 0.0_JPRB
        SP_OM(J,JLEV,YSP_OM%YVO%MP,IBL) = 0.0_JPRB
        SD_V3(J,JLEV,YSD_V3%YOTKE%MP,IBL) = 0.0_JPRB
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

CALL MPL_BARRIER(CDSTRING='SUINIF2')

!*       2.3   UPPER AIR
IF( NUMGPFLDS>0 .OR. LENCLD2 .OR. LEPCLD2 ) THEN
  CALL SUGRIDU(YDGEOMETRY,YDMODEL%YRML_PHY_SLIN%YREPHLI,YDMODEL%YRML_GCONF,YDMODEL%YRML_LBC,YDSP,YDGFL,KFILE)
ENDIF

CALL MPL_BARRIER(CDSTRING='SUINIF3')

!*       2.4   COMPUTE RT FROM T AND Q
!*       2.5   COMPUTE PD VD for NH case
IF (YDDYNA%LNHDYN) THEN
  CALL SPNH_CONV_NHVAR(YDGEOMETRY,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA,.FALSE.,YDSP)
ENDIF

!*       2.6   COMPUTE RT FROM T AND Q

IF (YDDYNA%LSPRT) THEN
  IF (LELAM) THEN
      CALL ESPECRT(YDGEOMETRY,YDGFL,YGFL,0,YDSP)
  ELSE
    IF (.NOT.(LIFSMIN .AND. LTRAJHR)) THEN
      CALL SPECRT(YDGEOMETRY,YDGFL,YGFL,0,YDSP)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    INITIALIZE CUMULATED FLUXES FROM INPUT FILE
!              -------------------------------------------

IF(LCUMFU.AND.LREACFU) THEN
  IF (LFBDAP) THEN
    CALL SUGRCFU(YDGEOMETRY,YDCFU,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,KFILE)
  ELSE
    IFILE=6
    CALL SUGRCFU(YDGEOMETRY,YDCFU,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,IFILE)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       4.    INITIALIZE INSTANTANEOUS FLUXES FROM INPUT FILE
!              -----------------------------------------------

IF(LXFU.AND.LREAXFU) THEN
  IF (LFBDAP) THEN
    CALL SUGRXFU(YDGEOMETRY,YDXFU,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,KFILE)
  ELSE
    IFILE=7
    CALL SUGRXFU(YDGEOMETRY,YDXFU,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,IFILE)
  ENDIF
ENDIF
CALL GSTATS(14,1)

IF (.NOT.(LIFSMIN .AND. LTRAJHR)) THEN
  CALL SPNORM(YDGEOMETRY,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA,YDSP,LDSPOR=LLSPOR)
ENDIF
CALL FLUSH(NULOUT)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUINIF',1,ZHOOK_HANDLE)
END SUBROUTINE SUINIF

