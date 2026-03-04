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

SUBROUTINE SU0YOMB(YDFPOS,YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDTCV_BGC,YDODB)

!**** *SU0YOMB*  - INITIALIZE LEVEL 0 COMMONS AND SOME HIGHER (PART 2)

!     PURPOSE.
!     --------
!           INITIALIZE LEVEL 0 COMMONS (CONSTANT ALONG ALL THE JOB)
!        AND SOME HIGHER LEVEL COMMONS. CALLS ROUTINES THAT PERFORM
!        ALL THE PREPARATIONS NEEDED TO EXECUTE MODEL. THE TASK OF
!        INITIALIZING THE COMMONS IS DIVIDED BETWEEN TWO ROUTINES
!        (SU0YOMA AND SU0YOMB) TO AVOID PROBLEMS WITH USING POINTER
!        ARRAYS WHOSE DIMENSIONS ARE NOT DEFINED UNTIL AFTER CALLING SUDIM.

!        ky: some constraints to be known for future code reorganisation (this list is not comprehensive):
!         - SUPPVI does not depend on geometry and could be called earlier in the set-up (somewhere in SU0YOMA)
!         - SUPP must be called after SUGFL and partly uses horizontal geometry.
!         - SUIOS does not depend on geometry and could be called earlier in the set-up (somewhere in SU0YOMA)
!         - SUPTRTC does not depend on geometry.
!         - SUDFI could probably be called earlier; requires to know TSTEP.
!         - SUSPNG could be called inside SUDYN, like calls to SU(E)HDF and SURAYFRIC (diffusive processes).
!         - SUVAREPS must be called before SUGRIB.

!     INTERFACE.
!     ----------
!        *CALL* *SU0YOMB*

!        EXPLICIT ARGUMENTS
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!      see below

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!        ORIGINAL      : 87-10-15

!     MODIFICATIONS.
!     --------------
!        R. El Khatib  : 31-Aug-2007 Initialise uninitialised variable
!        B. Strajnar   : 06-02-08 Allow call to SUINFCE for LSPFCE=.F.
!        B. Chapnik    : 15-06-08 Add setup of Jk modulation wrt/ level and parameter
!        R. El Khatib  : 24-Oct-2008 Merge sueoph and suoph
!        K. Yessad     : 09-09-08 No recalc of Gaussian lat and weights in SUDIL
!        K. Yessad     : 15-09-08 Prune conf 951.
!        A. Deckmyn    : 01-10-08 LAM wavelets
!        P. Termonia   : 22-12-08 move the call to SUDFI here (from SUINI)
!        A. Dethof     : 17-11-08 LOG NOX JB
!        K. Yessad (Aug 2009): add call to SUPPVI.
!        K. Yessad (Aug 2009): prune conf 912, externalise conf 911.
!        K.Yessad (Feb 2010): use YM_RADTC and RFORADTC
!        R. El Khatib 13-Jul-2010 Move surfex allocations below suphmse
!        R. El Khatib 13-Jul-2010 Move surfex allocations below suphmse
!        K. Yessad (Sep 2010): organigramme simplification.
!        K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!        P. Marguinaud 01-Jan-2011 IO server setup (send parameters with MPI)
!        G. Kerdraon   : Feb 2011 Call SU_GRIB_API if NOT LELAM
!        K. Yessad (jul 2011): reorder calculations in order to use new structures in LAM models.
!        T.Wilhelmsson (Aug 2011) SUSC2B => SUSC2B + SUSC2C
!        M. Fisher   7-March-2012 Move Jb setup here (from SU0YOMA)
!        K. Yessad (dec 2011): various contributions.
!        R. El Khatib 09-Mar-2012 : Unconditional call to sualdyn_ddh
!        R. El Khatib 23-Mar-2012 : Fix bounds checking issue
!        B. Bochenek (Apr 2012): call to SUNDDH moved from SU0YOMA
!        R. El Khatib 26-Jul-2012 : SUVV1 (previously part of SUVERT)
!                                 + SUFPG (previously part of SUBFPOS)
!        M.Fisher 15-Feb-2013 Introduce CVA_STRUCT, SCALP_STRUCT, CVA_DATA
!        K. Yessad (july 2013): various modifications for geometry set-up.
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        D. Degrauwe, 2013-11 : Setup of flexible interface (INTFLEX).
!        R. El Khatib  04-Jul-2014 call sualspa1 before sualspa to facilitate the
!         subsequent allocation of the spectral structure by alloc_spec
!        K. Yessad (july 2014): some reorganisation in set-up.
!        R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!        P. Marguinaud 10-Oct-2014 Add LGRIB_API
!        R. El Khatib : 03-Dec-2014 skeleton of the configuration 903
!        A. Geer      27-Jul-2015   VarBC is now an object, for OOPS
!        P. Lopez     Dec 2015   Skip Jb allocations for configuration 501 (TL test).
!        M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!        SJ Lock : Jan-2016  Cleaning SPPT routines
!        R. El Khatib 17-Aug-2016 move suoph up to su0yoma
!        O.Marsden Aug 2016  Removed use of SPA3
!        K. Yessad (Dec 2016): Prune obsolete options.
!        B.Bochenek(Feb 2017): Fix for 601 - CALL SUSPSDT before SUSCAL/SUESCAL
!        B.Bochenek(Feb 2017): Temporary fix for AROME in ALLOCATE_SPEC
!        J. Hawkes 22-Nov-2017 Initialize part of wave model for IO server
!        Y. Michel, MF, June 2018 Extension of the control variable for sqrt EnVar
!        S. Massart   19-Feb-2019 Augmented control variable
!        Y. Michel, MF,  Mar 2019 Extention of the control variable for sqrt EnVar
!        C. Lupu 29-Mar-2019 Allow call to SUINSKTE
!        M. Leutbecher Oct-2020 SPP abstraction
!     ------------------------------------------------------------------

USE TYPE_MODEL      , ONLY : MODEL
USE GEOMETRY_MOD    , ONLY : GEOMETRY
USE VARIABLES_MOD   , ONLY : VARIABLES, VARIABLES_CREATE, VARIABLES_DELETE
USE FIELDS_MOD      , ONLY : FIELDS, FIELDS_CREATE, FIELDS_DELETE, FIELDS_CONTAIN
USE MTRAJ_MOD       , ONLY : MTRAJ
USE PARKIND1        , ONLY : JPIM, JPRB
USE YOMHOOK         , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN          , ONLY : NULOUT
USE JO_TABLE_MOD    , ONLY : JO_TABLE
USE YOMCT0          , ONLY : LR2D ,NCONF, LBACKG, LSPBSBAL, LELAM,LARPEGEF,&
 &                           LOBS, LOBSC1, LSCREEN, LIFSTRAJ, LIFSMIN, NFPOS,&
 &                           NFRCO, LGRIB_API
USE YOMARG          , ONLY : NGRIBFILE
USE YOMVAR          , ONLY : LMODERR , LVARBC, LJCDFI, LMONITOR_FCDEPAR,LBACKGERENORM, LECV,LENSCV
USE ALGORITHM_STATE_MOD  , ONLY : GET_NUPTRA
USE YOMJG           , ONLY : JB_STRUCT
USE YOMJQ           , ONLY : YGERRMOD, LSTATMERR
USE YOMLCZ          , ONLY : L_SUBSPACE_SVS, LFORCEWR
USE YOMCVA          , ONLY : CVA_STRUCT, SCALP_STRUCT, CVA_DATA
USE YOMSKTER        , ONLY : YGSKTER
USE TRAJECTORY_MOD  , ONLY : LTRAJRESET, LTRAJHR, ALLOCATE_TRAJECTORY,&
 &                           READ_TRAJECTORY, GET_TRAJ_GRID, LREADGPTRAJ, LTRAJHR_ALTI
USE YOMGWDIAG       , ONLY : SETUP_GWDIAG
USE YEMJK           , ONLY : LEJK
USE YOMINI          , ONLY : LDFI
USE YOE_CUCONVCA    , ONLY : INI_CUCONVCA
USE MODULE_RADTC_MIX, ONLY : YM_RADTC ,SUPTRTC




USE YEMLBC_MODEL    , ONLY : SUELBC_INIT

USE YOMIO_SERV      , ONLY : IO_SERV_C001
USE YOM_GRIB_CODES  , ONLY : NGRBNOXLOG
USE SPNG_MOD        , ONLY : SUSPNG
USE YOMTRAJ         , ONLY : LTRAJALLOC, TRAJEC
USE YOMMP0          , ONLY : NPROC
USE TESTVAR_MIX     , ONLY : SETUP_TESTVAR
USE VARBC_CLASS     , ONLY : CLASS_VARBC
USE TOVSCV_MOD      , ONLY : TOVSCV
USE TOVSCV_BGC_MOD  , ONLY : TOVSCV_BGC
USE YOMSPJB         , ONLY : BACKGROUND, ALLOCATE_JB_REF_STATE
USE YOMJBECV        , ONLY : LJB_ALPHA_CV, LECPHYSPARECV, READ_FG_ECV
USE YOMJBPAR1DECV   , ONLY : SUPARECVMIN
USE YOMJBECPHYSECV   ,ONLY : SUINFCE_ECPHYS, LSOLARCST
USE MPL_MODULE      , ONLY : MPL_END, MPL_BARRIER
USE IOSTREAM_MIX    , ONLY : IOSTREAM_STATS,YGBH
USE DBASE_MOD   , ONLY : DBASE

USE FULLPOS         , ONLY : TFPOS
USE YOMFP_SERV, ONLY : FP_SERV_C001
USE CONTROL_VECTORS_MOD
USE SPECTRAL_FIELDS_MOD
USE YOMMODERR       , ONLY : SPCTLMODERR
USE OBS_STORE_OPTIONS_MOD, ONLY : YDOBS_STORE_OPTIONS
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TFPOS),       INTENT(OUT)   :: YDFPOS
TYPE(GEOMETRY),    INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS),      INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ),       INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL),       INTENT(INOUT) :: YDMODEL
TYPE(JO_TABLE),    INTENT(INOUT) :: YDJOT
TYPE(CLASS_VARBC), INTENT(INOUT) :: YDVARBC
TYPE(TOVSCV),      INTENT(INOUT) :: YDTCV
TYPE(TOVSCV_BGC),  INTENT(INOUT) :: YDTCV_BGC
CLASS(DBASE),      INTENT(OUT)   :: YDODB

CHARACTER (LEN = 35) ::  CLINE
CHARACTER(LEN=40) :: CLFILE
CHARACTER (LEN=3)  :: CLMAX

INTEGER(KIND=JPIM) :: ICONF

LOGICAL :: LLDIMO, LL_ALLOC_RLANBUF

INTEGER(KIND=JPIM) :: JGFL, ISTEP, JSTGLO, ICEND, IBL, ISIZEG, IWINLEN
INTEGER(KIND=JPIM) :: IGRIB(YDMODEL%YRML_GCONF%YRDIMF%NS3D+YDMODEL%YRML_GCONF%YRDIMF%NS2D)
LOGICAL :: LLASTRAJ
REAL(KIND=JPRB), ALLOCATABLE :: ZGMV5(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGMV5S(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGFL5(:,:,:,:)

TYPE(CONTROL_VECTOR) :: YLTEMP
TYPE(VARIABLES)      :: YL_VARS, YNLVARS
TYPE(FIELDS)         :: YL_TRAJ, YLINFCE
!!$CHARACTER(LEN=9)   :: CLCONF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

#include "get_spp_conf.intfb.h"
#include "gstats_output_ifs.intfb.h"
#include "ini_spp.intfb.h"
#include "sujbwavelet0.intfb.h"
#include "suecges.intfb.h"
#include "suejbbal.intfb.h"
#include "suejbcov.intfb.h"
#include "suejk.intfb.h"
#include "suejknorm.intfb.h"
#include "suelges.intfb.h"
#include "suelljk.intfb.h"
#include "suemodjk.intfb.h"
#include "suescal.intfb.h"
#include "suevargp.intfb.h"
#include "sueinfce.intfb.h"
#include "suensdim.intfb.h"
#include "suenscov.intfb.h"
#include "suejbwavelet.intfb.h"
#include "suejbwavelet_bmatrix.intfb.h"
#include "sualdyn_ddh.intfb.h"
#include "sualgco.intfb.h"
#include "sualmdh.intfb.h"
#include "suspvariables.intfb.h"
#include "sualspa.intfb.h"
#include "sualnud.intfb.h"
#include "sualtdh.intfb.h"
#include "subfpos.intfb.h"
#include "sucfu.intfb.h"
#include "sudyn.intfb.h"
#include "sugrib.intfb.h"
#include "su_grib_api.intfb.h"
#include "suiau.intfb.h"
#include "suinterpolator.intfb.h"
#include "suios.intfb.h"
#include "suiostream.intfb.h"
#include "sulsforc.intfb.h"
#include "sumcclag.intfb.h"
#include "sumddh.intfb.h"
#include "sunddh.intfb.h"
#include "suphy.intfb.h"
#include "supp.intfb.h"
#include "surand1.intfb.h"
#include "sures.intfb.h"
#include "surlx.intfb.h"
#include "susc2b.intfb.h"
#include "susc2c.intfb.h"
#include "susimpr.intfb.h"
#include "suspsdt.intfb.h"
#include "suvareps.intfb.h"
#include "suxfu.intfb.h"
#include "sumcuf.intfb.h"
#include "io_serv_suiosctmpl.intfb.h"
#include "allocate_empty_trajectory.intfb.h"
#include "suintflex.intfb.h"
#include "fp_serv_suiosctmpl.intfb.h"
#include "fp_serv_cpfpfilter.intfb.h"


!!$#include "transinvh.intfb.h"
!!$#include "suinif.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0YOMB',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL5=>YDMTRAJ%YRGFL5,YDGMV5=>YDMTRAJ%YRGMV5, YDGFL=>YDFIELDS%YRGFL,YDGMV=>YDFIELDS%YRGMV, &
 & YDSURF=>YDFIELDS%YRSURF, &
 & YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM,  YDMP=>YDGEOMETRY%YRMP, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDML_LBC=>YDMODEL%YRML_LBC, &
 & YDEWCOU=>YDMODEL%YREWCOU, &
 & YDRCOEF=>YDMODEL%YRML_PHY_RAD%YRRCOEF,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF, &
 & YDSPPT_CONFIG=>YDMODEL%YRML_GCONF%YRSPPT_CONFIG, YDSPPT=>YDMODEL%YRML_SPPT, &
 & YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG, YDSPP=>YDMODEL%YRML_SPP)

ASSOCIATE(NCHEM_ASSIM=>YGFL%NCHEM_ASSIM, NDIM5=>YGFL%NDIM5, YCOMP=>YGFL%YCOMP, &
 & NCMAX=>YDDIM%NCMAX, NGPBLKS=>YDDIM%NGPBLKS, NMSMAX=>YDDIM%NMSMAX, &
 & NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, NSMAX=>YDDIM%NSMAX, &
 & NUMP=>YDDIM%NUMP, &
 & NS2D=>YDDIMF%NS2D, NS3D=>YDDIMF%NS3D, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & NGPTOT=>YDGEM%NGPTOT, RSTRET=>YDGEM%RSTRET, &
 & YT5=>YDGMV5%YT5, &
 & MYMS=>YDLAP%MYMS, &
 & LSDDH=>YDLDDH%LSDDH, &
 & NALLMS=>YDMP%NALLMS, NPTRLL=>YDMP%NPTRLL, NPTRMS=>YDMP%NPTRMS, &
 & NUMLL=>YDMP%NUMLL, NPSURF=>YDMP%NPSURF, &
 & LRCOEF=>YDRCOEF%LRCOEF, &
 & LRAYSP=>YDSIMPHL%LRAYSP)
!     ------------------------------------------------------------------
CALL GSTATS(39,0)
CLINE='----------------------------------'

ICONF=NCONF/100

!*    Initialize post-processing
WRITE(NULOUT,*) '---- Set up post-processing ---------',CLINE
CALL SUPP(YGFL,YDMODEL%YRML_PHY_MF%YRPHY,YDMODEL%YRML_PHY_RAD%YREAERATM)

!*    Initialize I/O-scheme
WRITE(NULOUT,*) '---- Set up I/O scheme --------------',CLINE
CALL SUIOS

!*    Initialize Aladin Jk handling
IF (LELAM.AND.(ICONF == 1)) THEN
 WRITE(NULOUT,*) '---- Set up Aladin Jk handling -',CLINE
 CALL SUEJK
ELSE
  LEJK=.FALSE.
ENDIF

!*    Initialize DDH (Horizontal domains diagnostics)
WRITE(NULOUT,*) '------ Set up DDH diagnostics --------',CLINE
CALL SUNDDH(YDGEOMETRY,YDMODEL)

IF (LBACKG .OR. LOBSC1 .OR. NCONF == 401 .OR. NCONF == 501 .OR. NCONF == 601 .OR. NCONF == 801 .OR. NCONF == 701) THEN

  IF (.NOT.LR2D) THEN
    WRITE(NULOUT,*) '------ Set up Jb parameters ------------',CLINE
    ALLOCATE(JB_STRUCT)
    ALLOCATE(CVA_DATA)
 CALL SUJB(YDGEOMETRY,YDMODEL%YRML_GCONF%YRDIMECV,YGFL,JB_STRUCT,CVA_DATA)
  ENDIF
ENDIF

WRITE(NULOUT,*) ' LIFSMIN = ', LIFSMIN
IF(ASSOCIATED(JB_STRUCT)) THEN
  LLDIMO=.NOT. JB_STRUCT%WJBCONF%LJBWSTATS .OR. .NOT. LIFSMIN
  WRITE(NULOUT,*) ' JB_STRUCT%WJBCONF%LJBWAVSTATS = ', JB_STRUCT%WJBCONF%LJBWSTATS
ELSE
  LLDIMO=.NOT. LIFSMIN
ENDIF
!*     Allocate spectral arrays
WRITE(NULOUT,*) '-- Set up spectral arrays allocation-',CLINE
CALL SUSPVARIABLES(YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA%LNHX,&
 & YDMODEL%YRML_DYN%YRDYNA%LNHDYN)
IF(NCONF /= 901.AND.NCONF /= 923) THEN
  IGRIB(1:NS3D)=YDMODEL%YRML_GCONF%YRDIMF%NGRBSP3(:)
  IGRIB(NS3D+1:NS3D+NS2D)=YDMODEL%YRML_GCONF%YRDIMF%NGRBSP2(:)
  WRITE(NULOUT,*)'-- Calling ALLOCATE_SPEC ',IGRIB(1:NS3D+NS2D)
!!  CALL ALLOCATE_SPEC(YDFIELDS%YRSPEC, NFLEVL, NFLEVG, NUMP, MYMS, NSMAX, NMSMAX, NALLMS,&
  CALL CREATE_SPEC(YDFIELDS%YRSPEC, NFLEVL, NFLEVG, NUMP, MYMS, NSMAX, NMSMAX, NALLMS,&
                   & NPTRMS, NUMLL, NPTRLL, NPSURF, NS3D, NS2D, IGRIB)
  CALL SUALSPA(YDGEOMETRY,YDMODEL%YRML_DYN%YRDYNA%LGRADSP)
ELSE
  WRITE(NULOUT,*) 'SU0YOMB: No call for SUALSPA1'
ENDIF

!*    Allocate nudging arrays
WRITE(NULOUT,*) '-- Set up nudging arrays allocation-',CLINE
CALL SUALNUD(YDGEOMETRY,YGFL)

IF ( NCONF == 701 ) THEN





  CALL ABOR1('FC-Only/OIFS - CALL SUALCOS, SUALGES and SUELGES in su0yomb should never be called - Exit')

ENDIF

!*    Initialize Jb

SETUP_JB: IF (LBACKG .OR. LOBSC1 .OR. NCONF == 401 .OR. NCONF == 501 .OR. NCONF == 601 .OR. NCONF == 801) THEN

!     Setup cost function arrays
  WRITE(NULOUT,*) '------ Set up cost functions -----------',CLINE



  CALL ABOR1('FC-only/OIFS - CALL Sualcos in su0yomb should never be called - Exit')


!    Allocate Jb spectral arrays

  WRITE(NULOUT,*) '------ Allocate Jb spectral arrays -----',CLINE
  IF (NCONF==131 .OR. NCONF/100==4 .OR. NCONF/100==5 .OR. NCONF/100==6 .OR. NCONF/100==8 .OR.&
    & (LTRAJHR .AND. LIFSTRAJ)) THEN
    ALLOCATE(BACKGROUND)
  ENDIF
!!  CALL ALLOCATE_SPEC(JB_STRUCT%JB_DATA%SPJB, YDGEOMETRY, JB_STRUCT%CONFIG%SPVARS)
  CALL ALLOC_SPEC(JB_STRUCT%JB_DATA%SPJB, YDGEOMETRY, JB_STRUCT%CONFIG%SPVARS)
  IF(NCONF == 801)CALL SUALLR(YDGEOMETRY,JB_STRUCT)

!     Setup dimensioning for wavelet Jb
  WRITE(NULOUT,*) '------ Setup Wavelet Jb dimensioning ---',CLINE
IF (ASSOCIATED(JB_STRUCT)) THEN
  IF (JB_STRUCT%WJBCONF%LJBWAVELET) THEN
    IF (LELAM) THEN
      WRITE(NULOUT,*) '-- Setup Aladin wavelet Jb: preliminaries'
      WRITE(NULOUT,*) '(I have to do this before getting to the Control Vector setup!)'
      CALL SUEJBWAVELET(YDDIMV,JB_STRUCT)
      JB_STRUCT%WJBCONF%N_WAVELET_SCALES=0 ! We don't use it, so for safety, set it to 0
    ELSE
      WRITE(NULOUT,*) '-- Setup spectral filters for wavelet Jb'
      CALL SUJBWAVELET0(YDDIM,JB_STRUCT,CDFILE='wavelet.cv')
      WRITE(NULOUT,*) '-- Define multiple-resolution transforms for wavelet Jb'
      CALL SUJBWAVTRANS(YDDIM,JB_STRUCT)
    ENDIF
  ELSE
    JB_STRUCT%WJBCONF%N_WAVELET_SCALES=0
  ENDIF
ENDIF

  WRITE(NULOUT,*) '------- Allocate Jb arrays ----------',CLINE
  CALL SUALGES(YDGEOMETRY,JB_STRUCT,CVA_DATA)
  IF (LELAM) CALL SUELGES(YDGEOMETRY,JB_STRUCT)

!*    Initialize gridpoint buffers for analysis errors
  IF (NCONF /= 901) THEN
    WRITE(NULOUT,*) '---- Set up gridpoint buffers for analysis errors ---',CLINE
    CALL SUANEBUF(YDGEOMETRY,YGFL,JB_STRUCT)
  ENDIF
ENDIF SETUP_JB

!     Define pointers for transmission coefficients (simp. radiation).
WRITE(NULOUT,*) '-- Set up pointers for transmission coefficients-',CLINE
IF (LRCOEF) THEN
  CALL SUPTRTC(.TRUE.,YM_RADTC)
ELSE
  CALL SUPTRTC(.FALSE.,YM_RADTC)
ENDIF

!*    Initialize geometry parameters for gridpoint
!     error standard deviations in Aladin C+I zone (YEMVARGP)
IF (LELAM) THEN
  WRITE(NULOUT,*) '---- Set up Aladin geometry: YEMVARGP -',CLINE
  IF (ICONF == 1) CALL SUEVARGP(YDGEOMETRY)
ENDIF

!*    Set up DFI (we need geometry for SSDFI)
IF (LDFI .OR. LJCDFI) THEN
  WRITE(NULOUT,*) '---- Set up DFI: SUDFI -----',CLINE
  CALL SUDFI(YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_MF%YRPHY, &
    & YDMODEL%YRML_LBC%LTENC)
ENDIF

!*    Allocate OI CANARI grid points arrays
IF(ICONF == 7) THEN
  WRITE(NULOUT,*) '------- Allocate CANARI arrays ----',CLINE
  CALL SUALCAN(YDGEOMETRY)
ENDIF

!*    Allocate YOMGCO arrays.
IF(NFRCO /= 0) THEN
  WRITE(NULOUT,*) '---- Allocate YOMGCO arrays --------',CLINE
  CALL SUALGCO(YDGEM,YDMODEL%YRML_PHY_G%YRDPHY)
ENDIF

!*    Setup TESTVAR.
WRITE(NULOUT,*) '---- Set up TESTVAR ----------',CLINE
CALL SETUP_TESTVAR(YDMODEL%YRML_GCONF%YRRIP)

!     Initialize model:

!*    Initialize Dynamics
WRITE(NULOUT,*) '---- Set up model dynamics ----------',CLINE
CALL SUDYN(YDGEOMETRY,YDMODEL,NULOUT)

!*    Initialize vertical interpolator
WRITE(NULOUT,*) '---- Set up vertical interpolator -------',CLINE
CALL SUINTERPOLATOR(YDGEOMETRY,YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_DYN%YRSLINT)

!*    Initialize Relaxation
WRITE(NULOUT,*) '------ Set up Relaxation ',CLINE
CALL SURLX(YDDIM,YDDIMV,YDMODEL%YRML_GCONF%YRRIP,NULOUT)

!*    Initialize new sponge
WRITE(NULOUT,*) '---- Set up new sponge ----------',CLINE
CALL SUSPNG(YDMODEL%YRML_DYN%YRSPNG,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_DYN%YRDYNA,YDDIMV%NFLEVG,YDSTA%STZ)

!*    Initialize control of DFI: SUFW
IF (LDFI .OR. LJCDFI ) THEN
  WRITE(NULOUT,*) '---- Set up DFI initialization: SUFW ',CLINE
  CALL SUFW(YDLAP,YDDIM,YDGEOMETRY%YREGEO,YDGEOMETRY%YRELAP)
ENDIF

!*    Initialize Large Scale Forcings
WRITE(NULOUT,*) '---- Set up large scale forcings ----------',CLINE
CALL SULSFORC(YDMODEL%YRML_GCONF,NULOUT)


IF (.NOT.LR2D) THEN
  !*    Initialize Physics
  WRITE(NULOUT,*) '---- Set up model physics -----------',CLINE
  CALL SUPHY(YDGEOMETRY,YDMODEL,NULOUT)

! Fields for physics
  CALL YDFIELDS%YEC_PHYS_FIELDS%CREATE(YDGEOMETRY,YDMODEL%YRML_PHY_G%YRDPHY)

  !*    Initialize special keys for the climate version 2nd part
  WRITE(NULOUT,*) '---- Set up MCC climate model keys (lagged part) --',CLINE
  CALL SUMCCLAG(YDGEM,YDMODEL%YRML_GCONF,YDMODEL%YRML_AOC,YDMODEL%YRML_CHEM%YRCOMPO, &
       &  YDMODEL%YRML_CHEM%YRCHEM, YDMODEL%YRML_PHY_AER%YREAERSRC, &
       &  YDMODEL%YRML_PHY_EC%YREPHY, NULOUT, YDSURF=YDSURF)
ENDIF

!*    Initialize cumulated fluxes requests
WRITE(NULOUT,*) '------ Set up cumulated fluxes diags ---',CLINE
CALL SUCFU(YDGEOMETRY,YDFIELDS%YRCFU,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_RAD%YRERAD,YDMODEL%YRML_PHY_MF%YRPHY, &
 & NULOUT)

!*    Initialize instantaneous fluxes requests
WRITE(NULOUT,*) '------ Set up instantaneous fluxes diags ',CLINE
CALL SUXFU(YDGEOMETRY,YDFIELDS%YRXFU,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_MF%YRPHY,NULOUT)

!*    Initialize variables for VAREPS (NB: must be called before SUGRIB)
WRITE(NULOUT,*) '- Set up VAREPS configuration',CLINE
CALL SUVAREPS(YDMODEL%YRML_GCONF%YRRIP)

!*    Initialize filter for monitoring the coupling updates
WRITE(NULOUT,*) '------ Set up monitoring coupling-updates',CLINE
CALL SUMCUF(YDFIELDS%YMCUF,YDDIM,YDMODEL%YRML_GCONF%YRRIP)

!*    Initialise restart mechanism
WRITE(NULOUT,*) '---- Set up restart mechanism ------',CLINE
CALL SURES(YDMODEL%YRML_GCONF%YRRIP,NULOUT)

!*    Initialize buffers for gridpoint scanning, part B
IF (NCONF /= 901 .AND. NCONF /= 903) THEN
  WRITE(NULOUT,*) '---- Set up gridpoint scanning, part B ----',CLINE
  CALL SUSC2B(YDGEOMETRY,YDMODEL)
ENDIF

!*    Initialize forcing by coarser model: part A
WRITE(NULOUT,*) '--- Set up forcing by coarser model part A ---------',CLINE
CALL SUELBC_INIT(YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_LBC)

!*    Initialize buffers for gridpoint scanning, part C
IF (NCONF /= 901) THEN
  WRITE(NULOUT,*) '---- Set up gridpoint scanning, part C ----',CLINE
  CALL VARIABLES_CREATE(YNLVARS, .FALSE.)
  CALL SUSC2C(YDGEOMETRY,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_GCONF,YDMODEL%YRML_PHY_MF%YRPHY,&
   & YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_LBC%LTENC,YNLVARS,YDFIELDS%YRGFL, &
   & YDFIELDS%YRGMV,YDFIELDS%YRSURF)
  CALL VARIABLES_DELETE(YNLVARS)
ENDIF


!*    Initialize forcing by coarser model: part C
IF (LELAM) THEN
  WRITE(NULOUT,*) '--- Set up forcing by coarser model part C ---------',CLINE
  !OIFS - SUELBC_MODEL, _FIELDS, FIELDS_DIM exist in module
  !        yelbc_model.F90, maybe could dummy this mod but
  !        will have impact on other parts of code.
  CALL ABOR1('OIFS - LAM calls to SUELBC_MODEL and SUELBC_FIELDS in su0yomb.F90 not supported in OIFS - Exit')

ENDIF

!*    Initialize IAU handling
WRITE(NULOUT,*) '---- Set up IAU handling -',CLINE
CALL SUIAU(YDMODEL%YRML_GCONF%YRRIP)

!*    Initialize domains and masks for DDH

CALL SUALMDH(YDGEM,YDMODEL%YRML_DIAG)

IF ( LSDDH ) THEN
  WRITE(NULOUT,*) '---- Set up DDH diagnostic domains ',CLINE
  CALL SUMDDH(YDGEOMETRY,YDMODEL%YRML_DIAG)
ENDIF

IF (.NOT.LSDDH.AND.LRAYSP) THEN
  WRITE(NULOUT,*) '---- Set up for simp.rad.if not ddh',CLINE
  CALL SUSIMPR(YDGEOMETRY,YDMODEL%YRML_DIAG%YRMDDH)
ENDIF

!*    Initialize GRIB coding parameters
WRITE(NULOUT,*) '---- Set up files : GRIB parameters -',CLINE
CALL SUGRIB(YDDIM,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_PHY_G%YRDPHY,YDMODEL%YRML_PHY_MF%YRPHY)

!*    Full Post-processing (2nd part)
IF (NFPOS /= 0) THEN

  WRITE(NULOUT,*) '- Set up F-post processing, bundled part',CLINE
  CALL SUBFPOS(YDFPOS,YDGEOMETRY,YDMODEL,NFPOS)

  ! Objects to be subsequently constructed :

  !* Initialize GRIB API templates from Fullpos geometry
  IF (NGRIBFILE==1 .OR..NOT.LARPEGEF) THEN
    IF (LGRIB_API) THEN
      WRITE(NULOUT,'('' == Full-Pos constructor : setup GRIB API templates == '')')
      CALL SU_GRIB_API(YDGEOMETRY,YDFPOS%YFPVAB,YDMODEL%YRML_AOC%YRMCC%LMCC04,YGBH, &
                     & YDFPUSERGEO=YDFPOS%YFPGEOMETRY%YFPUSERGEO(1))
    ELSE
      WRITE(NULOUT,'(''Call to SU_GRIB_API switched off'')')
    ENDIF
  ENDIF

  !* Initialize the I/O server for Fullpos
  WRITE(NULOUT,'('' == Full-Pos constructor : setup io server file templates == '')')
ELSE
  !*    GRIB API
  IF (LGRIB_API) THEN
    WRITE(NULOUT,*) '- Set up GRIB API usage',CLINE
    CALL SU_GRIB_API(YDGEOMETRY,YDGEOMETRY%YRVAB,YDMODEL%YRML_AOC%YRMCC%LMCC04,YGBH)
  ELSE
    WRITE(NULOUT,*) 'Call to SU_GRIB_API switched off',CLINE
  ENDIF
ENDIF

! IOSTREAM
CALL SUIOSTREAM

!*    Memory allocation for cumulated DDH arrays (horizontal domains diags)
WRITE(NULOUT,*) '---- Set up DDH diags allocation --',CLINE
CALL SUALTDH(YDDIMV,YDMODEL%YRML_DIAG,YDMODEL%YRML_PHY_MF%YRARPHY,YDMODEL%YRML_PHY_MF%YRPHY)

!*    Memory allocation for dynamical DDH tendencies arrays
WRITE(NULOUT,*) '---- Set up dynamical DDH arrays allocation --',CLINE
CALL SUALDYN_DDH(YDGEOMETRY,YDMODEL%YRML_DIAG,YDMODEL%YRML_GCONF,&
 & YDMODEL%YRML_DYN%YRDYNA%LNHDYN,YDMODEL%YRML_DYN%YRDYNA%LNHX)

!     Setup model error arrays
IF (LMODERR.OR.LSTATMERR) THEN
  WRITE(NULOUT,*) '-- Setup model error arrays-',CLINE
  CALL SUMODERR(YDGEOMETRY,YDFIELDS%YRGMV,YDMODEL%YRML_GCONF)
ENDIF

!     Allocate space for control variable
IF (NCONF/100 /= 9 .AND. NCONF /= 1 .AND. NCONF /= 302 .AND. NCONF /= 201 .AND. NCONF /= 701) THEN
   IF(ASSOCIATED(JB_STRUCT)) THEN
     IF (.NOT.ASSOCIATED(CVA_DATA)) CALL ABOR1('Call SUALCTV, but CVA_DATA is not set up')
     IF (LECPHYSPARECV) THEN
       CVA_DATA%NVPARECV=YDMODEL%YRML_GCONF%YRDIMECV%NECV_1D
     ELSE
       CVA_DATA%NVPARECV=0
     ENDIF
     IF (LENSCV) THEN
       WRITE(NULOUT,*) '-- Allocate ens. control variable-',CLINE
       !* SUALCTV_ENS needs to know about ens. size in sqrt. Envar scheme
       CALL SUENSDIM(YDGEOMETRY)
       ALLOCATE(CTLVEC_STRUCT_ENS)
       !* change here if ensemble geometry is different
       CALL SUALCTV_ENS(YDGEOMETRY,CTLVEC_STRUCT_ENS,CVA_DATA)
     ENDIF
     WRITE(NULOUT,*) '-- Allocate static control variable-',CLINE
     ALLOCATE(CTLVEC_STRUCT)
     CALL SUALCTV(YDGEOMETRY,CTLVEC_STRUCT,CVA_DATA,JB_STRUCT,YDMODEL%YRML_GCONF%YRDIMECV)
   ELSE
     CALL ABOR1(' SU0YOMB: case where JB_STRUCT is used, but is not yet set up')
   ENDIF


ENDIF

!*    Initialize control of the Lanczos algorithm
LFORCEWR=.FALSE.
IF(.NOT.(ICONF == 0.OR.ICONF == 2.OR.NCONF == 302.OR.NCONF == 903)) THEN
  WRITE(NULOUT,*) '---- Set up Lanczos algorithm -------',CLINE
  LL_ALLOC_RLANBUF=(NCONF/100 == 6)
  IF(LL_ALLOC_RLANBUF) THEN
    CALL ALLOCATE_CTLVEC(YLTEMP)
    ISIZEG=YLTEMP%NSIZEG
    CALL DEALLOCATE_CTLVEC(YLTEMP)
  ELSE
    ISIZEG=0
  ENDIF
  CALL SULCZ(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,LL_ALLOC_RLANBUF,ISIZEG)
ENDIF

!*    Allocate the trajectory
IF(NCONF/100 == 1.OR.NCONF == 401.OR.NCONF == 501.OR.NCONF == 801.OR.NCONF == 601) THEN
  WRITE(NULOUT,*) '---- Allocate trajectory: SUALLT -------',CLINE
  CALL SUALLT(YDGEOMETRY,YDMTRAJ,YDFIELDS%YRGMV,YDFIELDS%YRSURF,YDMODEL)
  LTRAJRESET=.TRUE.
ENDIF

IF (LTRAJHR .AND. LIFSTRAJ) THEN
  WRITE(NULOUT,*) '---- Allocate trajectory: ALLOCATE_TRAJECTORY -------',CLINE
  CALL ALLOCATE_TRAJECTORY(YDGEOMETRY,YDGMV,YDGMV5,YDSURF,YDMODEL)
ENDIF

!*    Get trajectory values for NOX LOG variable
IF (LLDIMO) THEN !JEB avoid calculating this with wavelet
 IF(( NCHEM_ASSIM>0) .AND. NGRBNOXLOG> 0) THEN
  WRITE(NULOUT,*) '---- Get trajectory values for NOX LOG variable -------',CLINE
  IF (LIFSMIN.AND.(LTRAJHR.AND.LTRAJHR_ALTI) ) THEN

    ISTEP    = 0
    LLASTRAJ = .FALSE.
    CALL READ_TRAJECTORY(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,ISTEP,LLASTRAJ,LREADGPTRAJ)

    IF (LREADGPTRAJ) THEN

      ALLOCATE (ZGMV5(NPROMA,NFLEVG,YT5%NDIM,NGPBLKS))
      ALLOCATE (ZGMV5S(NPROMA,YT5%NDIMS,NGPBLKS))
      ALLOCATE (ZGFL5(NPROMA,NFLEVG,NDIM5,NGPBLKS))

      ZGMV5=0._JPRB
      ZGMV5S=0._JPRB
      ZGFL5=0._JPRB

      CALL GET_TRAJ_GRID(YDGEOMETRY,YDMODEL%YRML_GCONF,TRAJEC(0),YDGMV,YDGMV5,ZGMV5,ZGMV5S,ZGFL5,1)

      DO JGFL=1,NDIM5
        IF (YCOMP(JGFL)%LTRAJIO) THEN
          YDGFL5%GFL5(:,:,YCOMP(JGFL)%MP,:)=ZGFL5(:,:,YCOMP(JGFL)%MP,:)
          DO JSTGLO=1,NGPTOT,NPROMA
            ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
            IBL=(JSTGLO-1)/NPROMA+1
            YDGFL5%GFL5(1:ICEND,1:NFLEVG,YCOMP(JGFL)%MP,IBL)=ZGFL5(1:ICEND,1:NFLEVG,YCOMP(JGFL)%MP,IBL)
            YDGFL5%GFL5(ICEND+1:NPROMA,1:NFLEVG,YCOMP(JGFL)%MP,IBL)=0.0_JPRB
          ENDDO
        ENDIF
      ENDDO

      DEALLOCATE(ZGFL5)
      DEALLOCATE(ZGMV5)
      DEALLOCATE(ZGMV5S)

    ENDIF

  ENDIF
 ENDIF
ENDIF

!*    Allocate Jb linearisation state
IF (NCONF==131 .OR. NCONF/100==6 .OR. NCONF/100==8 .OR.&
  & (LTRAJHR .AND. LIFSTRAJ)) THEN
! Not needed for NCONF=401?
  WRITE(NULOUT,*) '---- Allocate Jb linearisation state -------',CLINE
  IF (.NOT.ASSOCIATED(BACKGROUND)) THEN
    CALL ABOR1 ('BACKGROUND HAS NOT BEEN ALLOCATED')
  ENDIF
  CALL ALLOCATE_JB_REF_STATE(YDGEOMETRY,YGFL,BACKGROUND,YDGMV,YDGMV5)
ENDIF

!*    Initialize Jb error covariance model
!*    Set up coefficients for humidity change of variable
IF (LBACKG .OR. LOBSC1) THEN
  WRITE(NULOUT,*) '---- Set up coefficients for humidity change of variable --',CLINE
  CALL SUJBCHVAR(YDVAB,YDDIMV,JB_STRUCT)
ENDIF

IF (LBACKG.OR.LSPBSBAL) THEN
  IF (LSPBSBAL.OR..NOT.JB_STRUCT%CONFIG%LJBENER) THEN
    WRITE(NULOUT,*) '---- Set up Jb balance operators --',CLINE
    IF(LELAM) THEN
      CALL SUEJBBAL(YDGEOMETRY,'STABAL96  ',JB_STRUCT)
    ELSE
      CALL SUJBBAL(YDGEOMETRY,'STABAL96  ',JB_STRUCT)
    ENDIF
  ENDIF
ENDIF

IF (LBACKG) THEN
  WRITE(NULOUT,*) '---- Set up Jb error covariances --',CLINE

!*   Copy trajectory into a FIELDS structure

  CALL VARIABLES_CREATE(YL_VARS, .TRUE.)
  CALL FIELDS_CREATE(YL_TRAJ,YDGEOMETRY,YDMODEL,YL_VARS)
  CALL VARIABLES_DELETE(YL_VARS)

  CALL GET_TRAJ_GRID(YDGEOMETRY,YDMODEL%YRML_GCONF,TRAJEC(0),YDGMV,YDGMV5,&
                   & YL_TRAJ%YRGMV%GMV,YL_TRAJ%YRGMV%GMVS,&
                   & YL_TRAJ%YRGFL%GFL,GET_NUPTRA())

  IF(LELAM) THEN
    IF (JB_STRUCT%WJBCONF%LJBWAVELET) THEN
      WRITE(NULOUT,*) 'Wavelet Jb in Aladin: good luck!'
      CALL SUEJBWAVELET_BMATRIX(YDDIMV)
    ELSE
      CALL SUEJBCOV(YDGEOMETRY,YDFIELDS,YL_TRAJ,JB_STRUCT)
      !* Envar : ens. data and localization setup
      IF (LENSCV) CALL SUENSCOV(YDGEOMETRY,YDFIELDS,YDMODEL)
    ENDIF
  ELSE
    IF (JB_STRUCT%WJBCONF%LJBWAVELET) THEN
      IF ((.NOT.JB_STRUCT%WJBCONF%LJBWSTATS).OR.(JB_STRUCT%WJBCONF%LHYBRID_JB)) THEN
        ! sujbwavelet must be called before sujbwavstats if LHYBRID_JB=T
        CLFILE='wavelet.cv'
        CALL SUJBWAVELET(YDGEOMETRY,YDFIELDS,YL_TRAJ,JB_STRUCT,CLFILE)
        CALL SUJBWAVELET_STDEVS(YDGEOMETRY,YDMODEL%YRML_GCONF,YDMODEL%YRML_CHEM%YRCHEM,JB_STRUCT)
      ENDIF

      IF (JB_STRUCT%WJBCONF%LJBWSTATS) THEN
        WRITE(NULOUT,*) '---- Calculate Wavelet Jb error covariances --',CLINE
        CALL SUJBWAVSTATS(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,JB_STRUCT)

       IF(LBACKGERENORM)THEN
       ! Option to reload the matrix and to compute renormalisation coeffs of wavelet B
          WRITE(CLMAX,'(I3)') NSMAX
          CLFILE='wavelet_out_T'//TRIM(ADJUSTL(CLMAX))//'.cv'
          CALL SUJBWAVELET(YDGEOMETRY,YDFIELDS,YL_TRAJ,JB_STRUCT,CLFILE)
          CALL SUJBWAVELET_STDEVS(YDGEOMETRY,YDMODEL%YRML_GCONF,YDMODEL%YRML_CHEM%YRCHEM,JB_STRUCT)
           WRITE(NULOUT,*) '---- Variational job: Compute renormalisation coefficient',&
         & 'for the variance -------------------'
          CALL SUJBWAVRENORM(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL%YRML_GCONF,JB_STRUCT)
       ENDIF

        ! The statistics file has been written and closed.
        ! Normal exit
        WRITE (NULOUT,*) 'Jb stats file has been written'
        WRITE (NULOUT,*) 'Normal exit'
        CALL IOSTREAM_STATS
        CALL GSTATS(0,1)
        CALL GSTATS_OUTPUT_IFS(YDMODEL%YRML_GCONF%YRRIP)
        IF (NPROC > 1) THEN
          CALL MPL_BARRIER(CDSTRING='SU0YOMB')
        ENDIF
        CALL MPL_END()
        STOP
      ENDIF

    ELSE
      IF (.NOT.ASSOCIATED(BACKGROUND)) THEN
        CALL ABOR1 ('BACKGROUND HAS NOT BEEN ALLOCATED')
      ENDIF
      CALL SUJBCOV(YDGEOMETRY,YDFIELDS,YDDIMF,YDMODEL%YRML_DYN%YRDYN,YL_TRAJ,JB_STRUCT)
      IF (NCONF==801.AND. .NOT.LBACKG) THEN
        CALL SUECGES(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,BACKGROUND,JB_STRUCT)
      ENDIF
    ENDIF
  ENDIF

  IF (LELAM) THEN
    IF (LEJK) THEN
      WRITE(NULOUT,*) '---- Set up Jk error covariances --',CLINE
      CALL SUELLJK(YDGEOMETRY,YDDIMF,NULOUT)
      CALL SUEJKNORM(YDDIMV,JB_STRUCT)
      WRITE(NULOUT,*) '---- Set up Jk modulation wrt level and param --',CLINE
      CALL SUEMODJK(YDSTA,YDDIMV,YDDIMF)
    ENDIF
  ENDIF

  IF(.NOT.(ICONF == 8) .AND..NOT.(ICONF == 6)) THEN
    WRITE(NULOUT,*) '---- Set up Background fields --',CLINE
    IF (.NOT.ASSOCIATED(BACKGROUND)) THEN
      CALL ABOR1 ('BACKGROUND HAS NOT BEEN ALLOCATED')
    ENDIF
    CALL SUECGES(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,BACKGROUND,JB_STRUCT)

    IF(.NOT.LELAM .OR. .NOT.JB_STRUCT%CONFIG%LSPFCE) THEN
      CALL VARIABLES_CREATE(YL_VARS, .FALSE.)
      CALL FIELDS_CREATE(YLINFCE,YDGEOMETRY,YDMODEL,YL_VARS)
      CALL VARIABLES_DELETE(YL_VARS)
!!$      IF(LTRAJHR) THEN
        CALL GET_TRAJ_GRID(YDGEOMETRY,YDMODEL%YRML_GCONF,BACKGROUND,YDGMV,YDGMV5,YLINFCE%YRGMV%GMV,&
         & YLINFCE%YRGMV%GMVS,YLINFCE%YRGFL%GFL,GET_NUPTRA())
!!$      ELSE
!!$        CALL SUINIF(YDGEOMETRY,YLINFCE%YRGFL,YLINFCE%YRSURF,0,YLINFCE%YRSPEC)
!!$        CLCONF(1:9)='0AAX00000'
!!$        CALL TRANSINVH(YLINFCE%GEOM,YLINFCE%YRGFL,YLINFCE%YRGMV,CLCONF,YLINFCE%YRSPEC)
!!$      ENDIF
        WRITE(NULOUT,*) '---- Set up Jb gridpoint background error stdev --',CLINE
        IF (LELAM) THEN
           CALL SUEINFCE(YDGEOMETRY,JB_STRUCT)
        ELSE
           CALL SUINFCE(YDGEOMETRY,YDMODEL%YRML_GCONF,YLINFCE,JB_STRUCT)
        ENDIF
      CALL FIELDS_DELETE(YLINFCE)
    ENDIF
    IF (.NOT.LELAM .AND. LECV .AND. IABS(NCONF/100)==1) THEN
      WRITE(NULOUT,*) '---- Set up JB extended control variable --',CLINE
      CALL READ_FG_ECV(YDGEOMETRY,YDMODEL%YRML_GCONF%YRDIMECV,YDFIELDS%YRSURF,LDBCK=.TRUE.)
      CALL READ_FG_ECV(YDGEOMETRY,YDMODEL%YRML_GCONF%YRDIMECV)
      IF (LECPHYSPARECV) THEN
        CALL SUINFCE_ECPHYS(YDGEOMETRY,JB_STRUCT)
        IF (LSOLARCST) CALL SUPARECVMIN(YDGEOMETRY)
      ENDIF
      IF (LJB_ALPHA_CV) CALL SETJBALPHACV(YDGEOMETRY, YDMTRAJ, YDFIELDS, YDMODEL, BACKGROUND, JB_STRUCT)
    ENDIF
  ENDIF

!*    Initialize Jq error covariances
  IF (LMODERR) THEN
    WRITE(NULOUT,*) '---- Set up Jq error covariances --',CLINE
    CALL SUJQ(YGERRMOD,YDGEOMETRY,SPCTLMODERR)
  ENDIF

!*    Delete temporary FIELDS-type trajectory structure
  CALL FIELDS_DELETE(YL_TRAJ)

ENDIF
!     Set up spectral stochastic diabatic tendencies
WRITE(NULOUT,'(A72)') '--- Set up stochastically perturbed parametrization tendencies '//CLINE
WRITE(NULOUT,*)       '      SPPT a.k.a. stochastic physics with spectral pattern'
CALL SUSPSDT(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YDSPPT_CONFIG,YDSPPT)


!*    Initialize scalar product
IF(.NOT.(NCONF == 1.OR.NCONF == 302.OR.ICONF == 2.OR.ICONF == 9.OR.ICONF == 7 )) THEN
  WRITE(NULOUT,*) '---- Set up scalar product --------',CLINE
  IF (LELAM) THEN
    ALLOCATE (SCALP_STRUCT)
    CALL SUESCAL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,SCALP_STRUCT,YDVARBC,JB_STRUCT,YDTCV)
  ELSE
    ALLOCATE (SCALP_STRUCT)
    CALL SUSCAL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,SCALP_STRUCT,YDVARBC,JB_STRUCT)
  ENDIF
ENDIF

!*    Initialize subspace
IF(.NOT.(ICONF == 0.OR.ICONF == 2.OR.NCONF == 302).AND.L_SUBSPACE_SVS) THEN
  WRITE(NULOUT,*) '---- Set up subspace for SV computation -------',CLINE
  CALL SU_SUBSPACE(YDGEOMETRY,YDDIMF)
ENDIF

!     Read in SKT forecast errors
IF (LSCREEN .AND. .NOT. LMONITOR_FCDEPAR) THEN
  WRITE(NULOUT,*) '---- Read in SKT EDA errors ---',CLINE
  CALL SUINSKTE(YGSKTER)
ENDIF

!     Read in forecast errors
IF (LSCREEN .AND. .NOT. LMONITOR_FCDEPAR) THEN
  WRITE(NULOUT,*) '---- Read in forecast errors ---',CLINE
  CALL INIFGER
ENDIF

!    If required allocate empty skeleton of trajectory structure
IF (.NOT. LTRAJALLOC .AND. ((NCONF /= 401).AND.(NCONF /= 501).AND. (&
                          & NCONF /= 601).AND.(NCONF /= 801) )) THEN
  CALL ALLOCATE_EMPTY_TRAJECTORY(YDDIM,YDMODEL%YRML_GCONF%YRRIP)
  ! ltrajalloc=.true.
ENDIF

WRITE(NULOUT,'(A72)') '--- Set up stochastic physics, SPBS, CABS '//CLINE
CALL INI_CUCONVCA(YDGEOMETRY,YDMODEL%YRML_DYN%YRDYNA,&
 & YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDMODEL%YRML_DYN%YRSL)
CALL SURAND1(YDGEOMETRY,YDMODEL%YRML_PHY_STOCH,YDMODEL%YRML_DYN%YRDYN,YDMODEL%YRML_GCONF%YRRIP, &
 &           YDMODEL%YRML_PHY_EC%YRECUCONVCA)
!     Set up stochastically perturbed parameterisation scheme

WRITE(NULOUT,*) '--- Set up stochastically perturbed parametrization scheme (SPP)  ',CLINE
CALL GET_SPP_CONF(YDMODEL%YRML_GCONF%YRRIP,YDSPP_CONFIG)
CALL INI_SPP(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YDSPP_CONFIG,YDSPP)

!     provide optional GEOMETRY argument here, since these are MODEL processes rather than IO_SERV ones
WRITE(NULOUT,*) '---- Create io server file templates --',CLINE
! Set up gravity wave diagnostics
CALL SETUP_GWDIAG(YDDIM)
!    Set up flexible physics-dynamics interface
WRITE(NULOUT,*) '---- Set up flexible physics-dynamics interface --',CLINE
CALL SUINTFLEX(YGFL,YDMODEL%YRML_PHY_MF%YRARPHY,YDMODEL%YRML_PHY_MF%YRPHY)
CALL FIELDS_CONTAIN(YDFIELDS,YDGEOMETRY,YDMODEL)

WRITE(NULOUT,*) '-------------------------------------',CLINE
WRITE(NULOUT,*) '------ END OF SETUPS at level 0 -----',CLINE
WRITE(NULOUT,*) '-------------------------------------',CLINE
WRITE(NULOUT,*) ' '
WRITE(NULOUT,*) ' '
CALL FLUSH(NULOUT)
CALL GSTATS(39,1)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SU0YOMB',1,ZHOOK_HANDLE)
END SUBROUTINE SU0YOMB
