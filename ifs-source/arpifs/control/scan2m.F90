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

#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE SCAN2M(YDGEOMETRY,YDFIELDS,YDGMV5,YDMODEL,CDCONF,LD_TST_GPGFL,LD_DFISTEP,YDSPEC,PTRAJEC,&
& YDGOM5,YDFPOS,YDFPDATA,YDODB)

!****-------------------------------------------------------------------
!**** *SCAN2M* - Grid-point space computations
!****-------------------------------------------------------------------
!     Purpose.   Computations in grid-point space
!     --------

!**   Interface.
!     ----------
!        *CALL* *SCAN2M (..)

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        CDCONF       : configuration of work (see doc.)
!        LD_TST_GPGFL : if T, do timestepping on grid-point GFL.

!     INPUT/OUTPUT:
!     -------------


!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud   *ECMWF*

! Modifications
! -------------
!   J. Vivoda  03-2002     PC schemes for NH dynamics (LPC_XXXX keys)
!   P. Smolikova 02-09-30 : variable d4 in NH (GPPAUX etc.)
!   G. Desroziers 02-09-16 : also treat Aladin mean wind trajectory
!   R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!   G.Mozdzynski: 02-10-01 support for radiation on-demand comms
!   C. Fischer   03-07-16: gppaux becomes part of gmv
!   O.Spaniel    03-04-15: cleaning-a same named entity from modules
!   R. El Khatib : 03-04-17 Fullpos improvemnts
!   C. Fischer : 03-06-03 - use igpcomp/ngptot_cap in call to cobs for LAM/Jo
!   R. El Khatib : 03-08-08 fullpos & gfl
!   G. Radnoti: 03-10-01 : update calls to eslextpol1
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   05-Jan-2004 Y. Seity   time stepping only for not coupled GFL(bug LGPQ=T)
!   M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!   27-Feb-2004 M. Tudor   time stepping only for during predictor
!   10-Jun-2004 J. Masek   NH cleaning (LPC_NOTR, LFULLIMP)
!   01-Jul-2004 K. Yessad  Make clearer the tests about PC schemes.
!   R. El Khatib : 04-10-21 Fullpos fix for extra gfl
!   Y. Seity : 04-11-16: new arguments to CTVTOT for R and Cp calc.
!   14-Mar-2005 R. Engelen   Fullpos of greenhouse gases
!   14-Mar-2005 A. Untch   Fullpos of aerosols
!   Y.Tremolet    01-Aug-2005 Number sections consistently with TL and AD
!   Y.Bouteloup : 05-01-19 : Add QL QI in the call to capotx
!   Y.Bouteloup : 05-01-31 : Add YR in the call to VPOS
!   G. Hello: 05-01-26 Fullpos fix for arome gfl
!   M.Hamrud  15-Jan-2005  Revised GPRCP
!   G. Desroziers and K. Yessad (sept 2005):
!    - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!    - adapt option LTRAJHR to METEO-FRANCE configurations.
!   R. El Ouaraini & R. El Khatib : 05-07-28 Allocate CUF fields and CALL fpinvtrcuf.
!   R. El Khatib  20-May-2005 NFPWIDE moved to YOMWFPDS
!   Y.Bouteloup : 05-12-25 : Add QR QS in the call to capotx
!   O.Spaniel :   22-02-2006 : Phasing CY31
!   M. Bellus :   29-Sep-2006 : Add prognostic convection variables (ALARO-0)
!    - YDAL, YDOM, YUAL, YUOM, YUEN, YUNEBH in the call to VPOS
!   E. Holm      05-10-20: Move jbvcoord_interpolate to cvargptl
!   08-08-2006  S. Serrar diagnostic tracers (TRAC) added to GFL fields
!   M.Hamrud      01-Jul-2006  Revised surface fields
!   K. Yessad (oct 2006): bug correction for LTRAJHR in part 2.
!   J.Haseler 27-Feb-2007  Generalise read of gridpoint trajectory
!   S. Serrar 22-03-2007: changes for the post-processing of ERA40 GFL fields
!   S. Serrar 17-07-2007 post-processing methane related GFL fields
!   S. Serrar 02-05-2008 test on LERA40 removed
!   Y. Seity  11-01-2008: bf for non-coupled GFL not initialized in E zone
!   R. El Khatib  17-Oct-2008 bf for ALADIN fullpos with equal distribution
!                             of C+I+E (RDISTR_E=1.)
!   K. Yessad Dec 2008: merge SLCOMM+SLCOMM1 -> SLCOMM.
!   K. Yessad Dec 2008: merge the different (E)SLEXTPOL.. -> (E)SLEXTPOL.
!   K. Yessad (Aug 2009): rewrite vertical interpolator in post-processing.
!   K. Yessad (Jan 2010): revised code (first draft)
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   J.Hague (Oct 2011): CALL COBSALL moved here from SCAN2MTL
!   J.Hague (Oct 2011): VSURF (Surface derived type), GFL, GMV, GMVS passed down
!   J.Hague (Oct 2011):-L_OK set to true after first iteration (moved from GP_MODEL_TL)
!   A. Geer  (Jan 2012): COBSALL before and after the GP model for diagnostic physics
!   A. Geer  (Mar 2012): GOM rewrite
!   G. Mozdzynski (May 2012): further cleaning
!   K. Yessad (Dec 2011)   : various contributions.
!   C. Soci  27-02-2012    : add Model orog and lsm for Mesan corr function
!   U. Andrae (Apr 2012)   : add graupels (ZGT) in the call to capotx
!   R. El Khatib 18-Jul-2012 Fullpos move away
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   F. Vana  28-Nov-2013   : Redesigned trajectory handling
!   R. El Khatib 02-Oct-2014 gp_model allocating large arrays on heap or stack
!   K. Yessad (July 2014): Move some variables.
!   R. El Khatib 12-Aug-2016 cleaning
!   E.Dutra/G.Arduini (Jan 2018) : Change in SP_SG, snow multi-layer
!   Y.Hirahara (Apr 2019) : Add SP_SL for CMEM
!   E.Holm (Jul 2019) : Change in interface to STORE_MAIN_TRAJ to allow writing
!                       out BG .ne. FG (needs YDGEOMETRY intent INOUT).
!-----------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : GPOPER, TSURF
USE YOMGMV             , ONLY : TGMV
USE YOMGFL             , ONLY : TGFL
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT, NULERR
USE YOMCT0             , ONLY : NCONF, LIFSMIN, LARPEGEF, &
 &                              LELAM, LIFSTRAJ, LOBSC1, NINTERPTRAJ,&
 &                              NOPT_MEMORY
USE YOMCT3             , ONLY : NSTEP
USE TRAJECTORY_MOD     , ONLY : LTRAJSAVE, LTRAJCST, LTRAJHR, LREADGPTRAJ,&
 &                              LTRAJHR_ALTI, LTRAJHR_SURF, STORE_MAIN_TRAJ
USE YOMLCZ             , ONLY : LFORCEWR, GPFORCEU, GPFORCEV, GPFORCET, GPFORCEQ, GPFORCESP
USE TESTVAR_MIX        , ONLY : LTESTINC
USE YOMECTAB           , ONLY : YECVAR
USE YOMTRAJ            , ONLY : NTRAJ_CST, MSTEPTRAJW, LPRTTRAJ, MTRAJ_GRIB,MTRAJ_GRIB2,MAIN_GRIB,&
 &                              NGP5, MIOTRAJSURF, MTYPE_SURF_TRAJ, LTRAJGP, TRAJ_TYPE
USE YOMOPH0            , ONLY : CFNTRAJHRSURF, LINC
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA
USE YOMSPJB            , ONLY : BACKGROUND
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE SUPERGOM_CLASS     , ONLY : CLASS_SUPERGOM
USE FULLPOS            , ONLY : TFPOS, TFPDATA
USE YOMPPC             , ONLY : NFPPHY,MFPPHY
USE YOMFPC             , ONLY : LFPPACKING, LOCEDELAY

USE DBASE_MOD, ONLY : DBASE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV5
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=9)    ,INTENT(IN)    :: CDCONF
LOGICAL             ,INTENT(IN)    :: LD_TST_GPGFL
LOGICAL             ,INTENT(IN)    :: LD_DFISTEP
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC
TYPE(TRAJ_TYPE)     ,INTENT(INOUT) :: PTRAJEC
TYPE(CLASS_SUPERGOM), OPTIONAL, INTENT(INOUT) :: YDGOM5
CLASS(DBASE)     ,INTENT(INOUT) :: YDODB
TYPE(TFPOS)         ,INTENT(IN), OPTIONAL    :: YDFPOS
TYPE(TFPDATA)       ,INTENT(INOUT), OPTIONAL :: YDFPDATA
!     ------------------------------------------------------------------
REAL(KIND=JPRB), POINTER :: ZGFL5(:,:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: ZTEMP(:,:)  ! for CST trajectory
!    GRIBEX cannot handle more than one type of REAL and thus cannot accept MKINDTRAJ arrays
REAL(KIND=JPRB), ALLOCATABLE :: ZTRAJ_BUF(:,:,:) ! *  intermediate trajectory buffer

REAL(KIND=JPRB),ALLOCATABLE :: ZPSP_SB(:,:,:,:), ZPSP_SG(:,:,:,:), ZPSP_SL(:,:,:),ZPSP_RR(:,:,:), ZPSP_CL(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZPSP_CI(:,:,:), ZPSP_X2(:,:,:), ZPSD_WS(:,:,:), ZPSD_VD(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZPSD_VX(:,:,:), ZPSD_VF(:,:,:), ZPSD_VV(:,:,:), ZPSD_VN(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZMCUFGP (:,:,:), ZCFUBUF (:,:,:), ZXFUBUF (:,:,:)
TYPE(TSURF) :: YLSURF
TYPE(TGMV)  :: YLGMV
TYPE(TGFL)  :: YLGFL

REAL(KIND=JPRB) :: ZDUMMY(YDGEOMETRY%YRDIM%NSPEC2)

INTEGER(KIND=JPIM) :: IBL, ICEND, IEND, IOFF, IST, JKGLO, JLEV, JL,&
 & IGPCOMP, JGFL, ISTEP, JSTEP, IDUMMY(1), IFPLAG

CHARACTER(LEN=20) :: CLFILE

LOGICAL :: LLTRAJHR, LLFIRST, LLAST, LLFPOS, LLOCEDELAY, LLOCE
LOGICAL :: LLCOBSALL, LL_GP_DONE, LL_UPD_DONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------


#include "abor1.intfb.h"



#include "ctvtot.intfb.h"
#include "gpnorm_gfl.intfb.h"
#include "gpnorm3.intfb.h"
#include "gp_model_stack.intfb.h"
#include "gp_model_heap.intfb.h"
#include "gp_model_heap2.intfb.h"
#include "gprcp.intfb.h"
#include "save_test4dinc.intfb.h"
#include "spnorm.intfb.h"
#include "write_grid_traj.intfb.h"
#include "gp_derivatives.intfb.h"
#include "fpwrncf.intfb.h"
#include "fullpos_drv.intfb.h"
#include "sp2gpmcuf.intfb.h"
#include "fullpos_precond.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SCAN2M',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL=>YDFIELDS%YRGFL,YDGMV=>YDFIELDS%YRGMV, YDSURF=>YDFIELDS%YRSURF, &
 & YDDIM=>YDGEOMETRY%YRDIM,  YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDCFU=>YDFIELDS%YRCFU, YDXFU=>YDFIELDS%YRXFU)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM5=>YGFL%NDIM5, NUMFLDS=>YGFL%NUMFLDS, &
 & YCOMP=>YGFL%YCOMP, YG=>YGFL%YG, YI=>YGFL%YI, YL=>YGFL%YL, YQ=>YGFL%YQ, &
 & YR=>YGFL%YR, YS=>YGFL%YS, YRSPEC=>YGFL%YRSPEC,&
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOT_CAP=>YDGEM%NGPTOT_CAP, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, &
 & GFL=>YDGFL%GFL, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, YT0=>YDGMV%YT0, &
 & LEGBRAD=>YDMODEL%YRML_PHY_EC%YREPHY%LEGBRAD, &
 & NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP, &
 & GFUBUF=>YDCFU%GFUBUF, XFUBUF=>YDXFU%XFUBUF, &
 & LAGPHY=>YDEPHY%LAGPHY, LEPHYS=>YDEPHY%LEPHYS, &
 & SD_VD=>YDSURF%SD_VD, SD_VF=>YDSURF%SD_VF, SD_VV=>YDSURF%SD_VV, SD_VX=>YDSURF%SD_VX, &
 & SD_WS=>YDSURF%SD_WS, SP_CI=>YDSURF%SP_CI, SP_CL=>YDSURF%SP_CL, SP_RR=>YDSURF%SP_RR, &
 & SP_SB=>YDSURF%SP_SB, SP_SG=>YDSURF%SP_SG, SP_SL=>YDSURF%SP_SL, SP_X2=>YDSURF%SP_X2, &
 & SD_VN=>YDSURF%SD_VN, &
 & YSD_VDD=>YDSURF%YSD_VDD, YSD_VF=>YDSURF%YSD_VF, YSD_VFD=>YDSURF%YSD_VFD, &
 & YSD_VVD=>YDSURF%YSD_VVD, YSD_VX=>YDSURF%YSD_VX, YSD_VXD=>YDSURF%YSD_VXD, &
 & YSD_WSD=>YDSURF%YSD_WSD, YSP_CI=>YDSURF%YSP_CI, YSP_CID=>YDSURF%YSP_CID, &
 & YSP_CLD=>YDSURF%YSP_CLD, YSP_RRD=>YDSURF%YSP_RRD, YSP_SBD=>YDSURF%YSP_SBD, &
 & YSP_SGD=>YDSURF%YSP_SGD, YSP_SLD=>YDSURF%YSP_SLD, YSP_X2D=>YDSURF%YSP_X2D, &
 & YSD_VND=>YDSURF%YSD_VND, YSD_VN=>YDSURF%YSD_VN, YDUPD=>YDSURF%YDUPD)

!     ------------------------------------------------------------------
!*       1.    INITIAL SETUP, MAIN LOOP CONTROL.
!     ------------------------------------------------------------------

!*       1.1   INITIAL SETUP

! Total number of grid points - ndglg*ndlon in global / ndguxg*ndlon in LAM
! used for model computations
IF (LELAM) THEN
  IGPCOMP=MIN(NGPTOT,NGPTOT_CAP)
ELSE
  IGPCOMP=NGPTOT
ENDIF

! Trajectory step
ISTEP =MSTEPTRAJW(NSTEP)

!     ------------------------------------------------------------------
!*       2.    INITIAL GRIDPOINT COMPUTATIONS
!     ------------------------------------------------------------------

IF (.NOT. LFORCEWR) CALL GSTATS(8,0)

IF (CDCONF(3:3) == 'C' .OR. CDCONF(3:3) == 'D') THEN
  CALL ABOR1('SCAN2M: CDCONF(3:3)= C or D not supported')
ENDIF

!J--- Z fields used in minimisation sby forcing (below) and by COBSALL (later)----
 IF ( CDCONF(6:6) == 'V'.OR.CDCONF(6:6) == 'C'.OR.CDCONF(6:6) == 'F' .OR. ((&
    & LTRAJHR.AND.LTRAJHR_SURF).AND.LTRAJSAVE.AND.LREADGPTRAJ) ) THEN
  IF(LIFSMIN) THEN
    ALLOCATE (ZGFL5(NPROMA,NFLEVG,NDIM5,NGPBLKS))
    IF (LTRAJGP) THEN
      CALL GSTATS(1051,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,ICEND,IBL)
      DO JKGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        IF(.NOT.ASSOCIATED(PTRAJEC%MAIN(IBL)%GFL))&
           & CALL ABOR1('READING TRAJEC%MAIN: TRAJEC%MAIN%GFL NOT ALLOCATED')
        ZGFL5(1:ICEND,:,1:NDIM5,IBL)=PTRAJEC%MAIN(IBL)%GFL(1:ICEND,:,1:NDIM5)
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1051,1)
      IF (LPRTTRAJ)  WRITE(NULOUT,*)'GREPTRAJ GET_TRAJ_GRID GRID istep=',ISTEP
    ELSE
      CALL ABOR1('TRAJEC%MAIN:NOT YET DONE/GFL')
    ENDIF
  ENDIF
ENDIF

! Forcing surface fields with interpolated trajectory
IF (LIFSMIN.AND.(LTRAJHR.AND.LTRAJHR_SURF).AND.LTRAJSAVE&
  & .AND. (ISTEP > 0)) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO,IEND,IBL,IST,ZTEMP)
  DO JKGLO=1,NGPTOT,NPROMA
    IST =1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    ALLOCATE(ZTEMP(NPROMA,NGP5))
    ZTEMP(IST:IEND,:) = PTRAJEC%SRFC(IBL)%MIX_SRFC(IST:IEND,:)
    CALL GPOPER(YDDIM,YDDYN,'SET0TOTRAJ',YDSURF,KBL=IBL,PFIELD=ZTEMP)
    DEALLOCATE(ZTEMP)
    IF (LPRTTRAJ.AND.PTRAJEC%SRFC(IBL)%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ FORCE TRAJEC%SRFC in SCAN2M'
  ENDDO
!$OMP END PARALLEL DO

ENDIF


! Forcing upper air fields with interpolated trajectory in the case LREADGPTRAJ=T.
IF (LIFSMIN.AND.(LTRAJHR.AND.LTRAJHR_ALTI).AND.LTRAJSAVE.AND.LREADGPTRAJ) THEN

  DO JGFL=1,NDIM5
    IF (YCOMP(JGFL)%LTRAJIO) THEN
      GFL(:,:,YCOMP(JGFL)%MP,:)=ZGFL5(:,:,YCOMP(JGFL)%MP,:)
    ENDIF
  ENDDO

  WRITE(NULOUT,*) 'NORMS IN SCAN2M AFTER HR TRAJECTORY FORCING'
  CALL SPNORM(YDGEOMETRY,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA,YDSPEC)
  CALL GPNORM_GFL(YDGEOMETRY,YDGFL)
ENDIF

!*       2.1  FORCE

IF (LFORCEWR) THEN
  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    DO JLEV=1,NFLEVG
      DO JL=1,ICEND
        GPFORCEU(JL,JLEV,IBL,0)=GMV(JL,JLEV,YT0%MU,IBL)
        GPFORCEV(JL,JLEV,IBL,0)=GMV(JL,JLEV,YT0%MV,IBL)
        GPFORCET(JL,JLEV,IBL,0)=GMV(JL,JLEV,YT0%MT,IBL)
        GPFORCEQ(JL,JLEV,IBL,0)=GFL(JL,JLEV,YQ%MP,IBL)
      ENDDO
    ENDDO
    DO JL=1,ICEND
      GPFORCESP(JL,IBL,0)    =GMVS(JL,YT0%MSP,IBL)
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('SCAN2M',1,ZHOOK_HANDLE)
  RETURN
ENDIF

!*       2.2 CONVERSION OF TTx ARRAY FROM 'TV' TO T

IF (YDDYNA%LSPRT) THEN
  CALL GSTATS(1033,0)
  IF (CDCONF(3:3)=='A'.OR.CDCONF(3:3)=='B'.OR.CDCONF(3:3)=='1') THEN
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO,ICEND,IBL)
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      CALL CTVTOT(YDGEOMETRY,YGFL,1,ICEND,GMV(1,1,YT0%MT,IBL),GFL(1,1,1,IBL))
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
  CALL GSTATS(1033,1)
ENDIF
IF(CDCONF(3:3) /= '0' .AND. YRSPEC%LGP) THEN
  DO JKGLO=1,NGPTOT,NPROMA
    IST =1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    CALL GPRCP(NPROMA,IST,IEND,NFLEVG,PGFL=GFL(1,1,1,IBL),&
     & PR=GFL(1,1,YRSPEC%MP,IBL))
  ENDDO
  CALL GPNORM3(YDGEOMETRY,GFL,NDIM,YRSPEC%MP,YRSPEC%CNAME)
ENDIF
IF (CDCONF(3:3) /= '0') THEN
  CALL GP_DERIVATIVES(YDGEOMETRY,YDFIELDS,YGFL)
ENDIF
!*       2.4 Save nonlinear fields for inner/outer loop tests
IF (LTESTINC .AND. NCONF==1 .AND. NSTOP>1 .AND. CDCONF(5:5)=='0') THEN
  CALL SAVE_TEST4DINC(YDGEOMETRY,YDGFL,YDGMV,YDGMV5,YDMODEL%YRML_GCONF,NSTEP)
ENDIF


!     ------------------------------------------------------------------
!*       3.    TRAJECTORY HANDLING
!     ------------------------------------------------------------------

! Store surface and grid-point trajectory
! ---------------------------------------

IF (LTRAJSAVE.AND.LTRAJCST) THEN
  ! This is a dirty fix.
  ! The trajectory should not contain any constants.
  IF (ISTEP==1) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IBL,IST,IEND,ZTEMP)
    DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      ! Temporary quantity ZTEMP required as the KIND of TRAJ_CST generaly differs
      !    from the one anticipated by GPOPER.
      ALLOCATE( ZTEMP(NPROMA,NTRAJ_CST))
      CALL GPOPER(YDDIM,YDDYN,'TRAJSTORECST',YDSURF,KBL=IBL,PFIELD=ZTEMP)
      PTRAJEC%CST(IBL)%MIX_CST(IST:IEND,:)=ZTEMP(IST:IEND,:)
!$OMP CRITICAL
      DEALLOCATE(ZTEMP)
!$OMP END CRITICAL
      IF (LPRTTRAJ.AND.PTRAJEC%CST(IBL)%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ TRAJEC%CST: constants saved in SCAN2M '
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
ENDIF


IF (NCURRENT_ITER == 0) THEN

  ! Saving main trajectory

  LLTRAJHR=LTRAJHR.AND.LTRAJHR_ALTI

  IF (LTRAJSAVE.AND. ((&
     & LLTRAJHR.AND.LIFSTRAJ).OR.(.NOT.LLTRAJHR.AND.LIFSMIN))) THEN

    ! Write out main upper air trajectory at inner-loop resolution
    ! in either spectral space (VOR,DIV,T,Q,(O3)) or in grid point
    ! space (Q,(O3)). Controlled by LREADGPTRAJ
    CALL STORE_MAIN_TRAJ(YDGEOMETRY,YDGMV,YDGMV5,YDMODEL,YDSPEC,GMV,GMVS,GFL,NSTEP,LREADGPTRAJ,PTRAJEC,BACKGROUND)

  ENDIF

  ! Saving surface trajectory

  LLTRAJHR=LTRAJHR.AND.LTRAJHR_SURF

  IF (LTRAJSAVE.AND. ((&
     & LLTRAJHR.AND.LIFSTRAJ).OR.(.NOT.LLTRAJHR.AND.LIFSMIN))) THEN

    ! Write out surface fields. This is always done in grid point space.
    IF (ISTEP > 0) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IBL,IST,IEND,ZTEMP)
      DO JKGLO=1,NGPTOT,NPROMA
        IBL=(JKGLO-1)/NPROMA+1
        IST =1
        IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        ALLOCATE(ZTEMP(NPROMA,NGP5))
        CALL GPOPER(YDDIM,YDDYN,'TRAJSTORE',YDSURF,KBL=IBL,PFIELD=ZTEMP)
        PTRAJEC%SRFC(IBL)%MIX_SRFC(IST:IEND,:)= ZTEMP(IST:IEND,:)
        DEALLOCATE(ZTEMP)
      ENDDO
!$OMP END PARALLEL DO

      IF (LIFSTRAJ) THEN
        IF (ISTEP == 1) THEN
          LLFIRST = .TRUE.
        ELSE
          LLFIRST = .FALSE.
        ENDIF
        LLAST=.TRUE.
        DO JSTEP=NSTEP+1,NSTOP
          IF(MSTEPTRAJW(JSTEP) > 0) LLAST=.FALSE.
        ENDDO

        IF(MIOTRAJSURF(ISTEP)== -1) THEN
          DO JSTEP=1,ISTEP-1
            IF(MIOTRAJSURF(JSTEP) /= JSTEP)THEN
              WRITE(NULERR,*)'HOLE IN SURF TRAJ',JSTEP,ISTEP,MIOTRAJSURF(JSTEP)
              CALL ABOR1('STORE_TRAJ_SURF:HOLE')
            ENDIF
          ENDDO

          CLFILE(1:20)=CFNTRAJHRSURF
          WRITE(CLFILE( 7: 8),'(I2.2)') GET_NUPTRA()
          WRITE(CLFILE(18:20),'(I3.3)') 0
          WRITE(NULOUT,*)'GREPTRAJ Saving surface trajectory TRAJEC%SRFC in ',CLFILE

          ALLOCATE(ZTRAJ_BUF(NPROMA,NGP5,NGPBLKS))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IBL)
          DO IBL=1,NGPBLKS
            ZTRAJ_BUF(:,:,IBL)=PTRAJEC%SRFC(IBL)%MIX_SRFC(:,:)
          ENDDO
!$OMP END PARALLEL DO

          !  Better to be rewritten into grib_api compliant way...
          CALL WRITE_GRID_TRAJ(YDGEOMETRY,YDRIP,CLFILE,MTRAJ_GRIB,IDUMMY,MTYPE_SURF_TRAJ,&
           & 0,NGP5,ZTRAJ_BUF,ZDUMMY,NINTERPTRAJ,LLFIRST,LLAST,MTRAJ_GRIB2)
          MIOTRAJSURF(ISTEP)=ISTEP
          DEALLOCATE(ZTRAJ_BUF)
        ELSE
          WRITE(NULOUT,*)'GREPTRAJ Warning: Surf traj at step,',NSTEP,ISTEP,&
           & 'already written'
        ENDIF

      ENDIF
      IF (LPRTTRAJ) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_SFC istep=',ISTEP,' in SCAN2M'
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!*       4.    GRIDPOINT MODEL.
!     ------------------------------------------------------------------

LLCOBSALL=.FALSE.
LL_GP_DONE = .FALSE.
IF (CDCONF(6:6) == 'V'.OR.CDCONF(6:6) == 'C'.OR.CDCONF(6:6) == 'F') THEN
  IF(LOBSC1) LLCOBSALL=.TRUE.
ENDIF
IF (NCONF==701) LLCOBSALL=.TRUE.

! GBRAD: Save SD_VN(YACCPR%MP) in SD_VN(YACCPR5%MP) to pass to GOM
IF(LEGBRAD)SD_VN(:,YSD_VN%YACCPR5%MP,:)=SD_VN(:,YSD_VN%YACCPR%MP,:)

IF (LLCOBSALL) THEN

  ! Take a copy of the surface fields, as they will be updated by GP_MODEL
  ! This is the first hacky step in creating a single call to cobsall
  ALLOCATE(ZPSP_SB(NPROMA,YSP_SBD%NLEVS,YSP_SBD%NDIM,NGPBLKS))
  ALLOCATE(ZPSP_SG(NPROMA,YSP_SGD%NLEVS,YSP_SGD%NDIM,NGPBLKS))
  ALLOCATE(ZPSP_SL(NPROMA,YSP_SLD%NDIM,NGPBLKS))
  ALLOCATE(ZPSP_RR(NPROMA,YSP_RRD%NDIM,NGPBLKS))
  ALLOCATE(ZPSP_CL(NPROMA,YSP_CLD%NDIM,NGPBLKS))
  ALLOCATE(ZPSP_CI(NPROMA,YSP_CID%NDIM,NGPBLKS))
  ALLOCATE(ZPSP_X2(NPROMA,YSP_X2D%NDIM,NGPBLKS))
  ALLOCATE(ZPSD_WS(NPROMA,YSD_WSD%NDIM,NGPBLKS))
  ALLOCATE(ZPSD_VD(NPROMA,YSD_VDD%NDIM,NGPBLKS))
  ALLOCATE(ZPSD_VX(NPROMA,YSD_VXD%NDIM,NGPBLKS))
  ALLOCATE(ZPSD_VF(NPROMA,YSD_VFD%NDIM,NGPBLKS))
  ALLOCATE(ZPSD_VV(NPROMA,YSD_VVD%NDIM,NGPBLKS))
  ALLOCATE(ZPSD_VN(NPROMA,YSD_VND%NDIM,NGPBLKS))

  ZPSP_SB=SP_SB
  ZPSP_SG=SP_SG
  ZPSP_SL=SP_SL
  ZPSP_RR=SP_RR
  ZPSP_CL=SP_CL
  ZPSP_CI=SP_CI
  ZPSP_X2=SP_X2
  ZPSD_WS=SD_WS
  ZPSD_VD=SD_VD
  ZPSD_VX=SD_VX
  ZPSD_VF=SD_VF
  ZPSD_VV=SD_VV
  ZPSD_VN=SD_VN

  LL_UPD_DONE = YDUPD%L_OK

ENDIF

!     ------------------------------------------------------------------
!*       5.    GRIDPOINT MODEL.
!     ------------------------------------------------------------------

IF (.NOT. LFORCEWR) CALL GSTATS(8,1)

! Post-processing first

LLFPOS=PRESENT(YDFPOS).AND.PRESENT(YDFPDATA)

IF (LLFPOS) THEN
  CALL GSTATS(30,0)
  IF (LAGPHY .AND. (NSTOP > 0) .AND. LEPHYS) THEN
    IFPLAG=1
  ELSE
    IFPLAG=0
  ENDIF
  LLOCEDELAY=LOCEDELAY .AND. (NSTOP > 0)
  LLOCE=.FALSE.
  ! monitoring of coupling-update frequency in gridpoint
  ALLOCATE(ZMCUFGP(YDGEOMETRY%YRDIM%NPROMA,YDFIELDS%YMCUF%NCUFNR,YDGEOMETRY%YRDIM%NGPBLKS))
  CALL SP2GPMCUF(YDGEOMETRY,YDFIELDS%YMCUF,ZMCUFGP)
  ! Post-processing itself
  IF (NSTOP > 0 .AND. LFPPACKING) THEN
    ! Prepare input data for post-processing, incl physical fields pre-packing
    CALL FULLPOS_PRECOND(YDGEOMETRY,YDFIELDS,YDMODEL,YLGMV,YLGFL,NSTEP,NSTOP,YLSURF,ZCFUBUF,ZXFUBUF)
    CALL FULLPOS_DRV(YDFPOS,YDGEOMETRY,YLGMV,YLGFL,YLSURF,YDFIELDS%YRCFU,YDFIELDS%YRXFU, &
     & ZCFUBUF,ZXFUBUF,YDFIELDS%YMCUF%NCUFNR,ZMCUFGP,YDMODEL,TSTEP,NSTEP,NSTOP,IFPLAG,LLOCEDELAY,LLOCE,YDFPDATA,NFPPHY,MFPPHY)
  ELSE
    CALL FULLPOS_DRV(YDFPOS,YDGEOMETRY,YDFIELDS%YRGMV,YDFIELDS%YRGFL,YDFIELDS%YRSURF,YDFIELDS%YRCFU,YDFIELDS%YRXFU, &
     & GFUBUF,XFUBUF,YDFIELDS%YMCUF%NCUFNR,ZMCUFGP,YDMODEL,TSTEP,NSTEP,NSTOP,IFPLAG,LLOCEDELAY,LLOCE,YDFPDATA,NFPPHY,MFPPHY)
  ENDIF
  DEALLOCATE(ZMCUFGP)
  IF (YDFPOS%YFPIOH%YNAMFPIOS%NFPWRITE==1) THEN
    IF (NSTOP == 0 .OR. .NOT.(LEPHYS.AND.LAGPHY)) THEN
      IF (LARPEGEF) THEN
        IF (LINC) THEN
          CALL FPWRNCF(YDFPOS%YFPCNT%CFPNCF,KSTEP=NSTEP,PTSTEP=TSTEP)
        ELSE
          CALL FPWRNCF(YDFPOS%YFPCNT%CFPNCF,KSTEP=NSTEP)
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  CALL GSTATS(30,1)
ENDIF

IF (.NOT. LFORCEWR) CALL GSTATS(8,0)

! Then the forecasting model

IF (CDCONF(4:4) /= '0') THEN
  IF (NOPT_MEMORY == 0) THEN
    CALL GP_MODEL_HEAP(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,PTRAJEC=PTRAJEC)
  ELSEIF (NOPT_MEMORY == 1) THEN
    CALL GP_MODEL_STACK(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,PTRAJEC=PTRAJEC)
  ELSE
    CALL GP_MODEL_HEAP2(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,PTRAJEC=PTRAJEC)
  ENDIF
  LL_GP_DONE = .TRUE.
ENDIF

IF (LLCOBSALL) THEN
  !!update YDPHYSMWAVE status as appropriate, to be used below COBSALL
  YDFIELDS%YEC_PHYS_FIELDS%YRPHYSMWAVE%LPHYS_MWAVE_FILLED_IN = .TRUE.

  CALL YDGOM5%MODEL_IN(&
    & YDGEOMETRY,YDGMV,YDSURF,YDMODEL%YRML_PHY_MF%YRPHY,YDMODEL%YRML_GCONF,&
    & YDMODEL%YRML_DYN%YRSLINT,YDMODEL%YRML_DYN%YRDYNA,&
    & NDIMGMV,YECVAR,YGFL,&
    & YDFIELDS%YEC_PHYS_FIELDS%YRPHYSMWAVE,GFL,GMV,GMVS,&
    & ZPSP_SB,ZPSP_SG,ZPSP_SL,ZPSP_RR,ZPSP_CL,ZPSP_X2,ZPSP_CI,&
    & ZPSD_VF,ZPSD_VV,ZPSD_VD,ZPSD_WS,ZPSD_VX,ZPSD_VN,NSTEP,&
    & LL_GP_DONE,LL_UPD_DONE)
  DEALLOCATE(ZPSP_SB,ZPSP_SG,ZPSP_SL,ZPSP_RR,ZPSP_CL,ZPSP_CI,ZPSP_X2,ZPSD_WS,ZPSD_VD,ZPSD_VX,ZPSD_VF,ZPSD_VV,ZPSD_VN)
ENDIF

 IF ( CDCONF(6:6) == 'V'.OR.CDCONF(6:6) == 'C'.OR.CDCONF(6:6) == 'F' .OR. ((&
    & LTRAJHR.AND.LTRAJHR_SURF).AND.LTRAJSAVE.AND.LREADGPTRAJ) ) THEN
  IF(LIFSMIN) THEN
    DEALLOCATE(ZGFL5)
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!*       7.    GRID POINT CALCULATIONS FOR ANALYSIS.
!     ------------------------------------------------------------------

!*       7.2  GRID POINT CALCULATIONS FOR CANARI   (cf. CA_SCAN2M)

IF (CDCONF(6:6) == '1') THEN



  ! OIFS - remove ca_scan2m since compile fails with gnu compiler and OIFS/FC-only
  !        build error "Explicit interface required for polymorphic argument at (1)"
  !        relates to use of types in the call
  CALL ABOR1('OIFS - ca_scan2m should not be called')

ENDIF

IF (.NOT. LFORCEWR) CALL GSTATS(8,1)

!     ------------------------------------------------------------------
!*       8.    "TIME-STEPPING" FOR GRID-POINT GFL FIELDS
!     ------------------------------------------------------------------

IF (LD_TST_GPGFL) THEN
  CALL GSTATS(1032,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,ICEND,IBL,JGFL,JLEV,JL)
  DO JKGLO=1,IGPCOMP,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    ICEND=MIN(NPROMA,IGPCOMP-JKGLO+1)
    DO JGFL=1,NUMFLDS
      IF(YGFL%YCOMP(JGFL)%LGP.AND.YGFL%YCOMP(JGFL)%LT1) THEN
        IF (YGFL%YCOMP(JGFL)%NCOUPLING==0) THEN
          YDFIELDS%YRGFL%GFL(1:ICEND,1:NFLEVG,YGFL%YCOMP(JGFL)%MP,IBL) =&
           & YDFIELDS%YRGFL%GFLT1(1:ICEND,1:NFLEVG,YGFL%YCOMP(JGFL)%MP1,IBL)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1032,1)
ENDIF
!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SCAN2M',1,ZHOOK_HANDLE)
END SUBROUTINE SCAN2M
