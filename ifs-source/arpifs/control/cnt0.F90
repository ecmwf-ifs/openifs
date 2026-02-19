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

SUBROUTINE CNT0(LDCOUPACTIVE,KCOMM)

!**** *CNT0*  - Routine which controls the job at level 0.

!     Purpose.
!     --------
!           Controls the job at level 0, the lowest level. The parameter
!       NCONF decides which path to follow.

!**   Interface.
!     ----------
!        *CALL* *CNT0

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M. Hamrud    : 09-08-01  Option removal
!      R. El Khatib : 01-08-07 Pruning options
!      F. Taillefer : 03-04-15 new config 932
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2005 Generalized IO scheme
!      M. Drusch    : 17-01-07  new config 302 / simplified EKF for soil moisture
!      K. Yessad    : 15-Sep-08 prune conf 951.
!      F. Taillefer : 10-04-01 change config 931 (main calling routine)
!      K. Yessad    : 27-Apr-10 prune conf 903 and 940.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib : 08-Jul-2014 Model fields should not be deallocated before the IO server has completed its tasks.
!      K. Yessad (July 2014) Encapsulate YOMRIP.
!     F. Vana  05-Mar-2015  Support for single precision
!      B. Bochenek (Apr 2015): Phasing: add call to geometry_unset
!      R. El Khatib : 03-Dec-2014 skeleton of the configuration 903
!      R. El Khatib : 01-Dec-2014 EC_DATE_AND_TIME + memory statistics control
!      A. Geer      27-Jul-2015   VarBC is now an object, for OOPS
!      R. El Khatib 10-Dec-2015 NSTEP  argument to opdis
!      P. Lean      17-Aug-2016   ODB is now an object, for OOPS
!      P. Lean      22-Mar-2017   Jo-table is now an object, for OOPS
!      R. El Khatib  04-Jun-2018 refactor suct1 against monio
!      R. El Khatib  03-Sep-2018 new configuration 904 which is a test program for change of resolution of an object FIELDS
!      P. Lopez     12-Oct-2018   Passed YDFPOS to CTL1
!       2019-04-26, Adrien Napoly and J.M. Piriou: 933 configuration operates both "old" 931 & 932.
!     ------------------------------------------------------------------

USE PARKIND1        , ONLY : JPRD, JPIM, JPRB, JPIB
USE YOMHOOK         , ONLY : LHOOK, DR_HOOK, JPHOOK
USE JO_TABLE_MOD    , ONLY : JO_TABLE
USE YOMCT0          , ONLY : NCONF, LECMWF, LELAM
USE YOMCT3          , ONLY : NSTEP
USE YOMLUN          , ONLY : NULOUT
USE YOMMP0          , ONLY : NOUTTYPE, LSLDEBUG, NPROC, LOUTPUT, LMPOFF, MYPROC
USE YOMVAR          , ONLY : LFDBERR
USE MPL_MODULE      , ONLY : MPL_BARRIER
USE YOMSIG          , ONLY : RESTART
USE EINT_MOD        , ONLY : UNUSED_HALO_STATS
USE MODEL_MOD       , ONLY : MODEL, MODEL_SET, MODEL_DELETE
USE FIELDS_MOD      , ONLY : FIELDS, FIELDS_DELETE
USE MTRAJ_MOD       , ONLY : MTRAJ, MTRAJ_DELETE
USE GEOMETRY_MOD    , ONLY : GEOMETRY, GEOMETRY_DELETE
USE IOSTREAM_MIX    , ONLY : CLOSE_IOSTREAM, Y_IOSTREAM_FDB, IOSTREAM_STATS
USE VARBC_CLASS     , ONLY : CLASS_VARBC
USE TOVSCV_BGC_MOD  , ONLY : TOVSCV_BGC
USE TOVSCV_MOD      , ONLY : TOVSCV
USE DBASE_MOD       , ONLY : DBASE



USE FULLPOS         , ONLY : TFPOS
#ifdef WITH_ATLAS
USE ATLAS_MODULE       , ONLY : ATLAS_LIBRARY
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT(IN), OPTIONAL :: LDCOUPACTIVE
INTEGER(KIND=JPIM),  OPTIONAL :: KCOMM

INTEGER(KIND=JPIM) ::  ICONFI
INTEGER(KIND=JPIM) ::  ITER

LOGICAL :: LLFIRSTCALL=.TRUE.

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE(GEOMETRY)     :: YRGEOMETRY
TYPE(MODEL),TARGET :: YRMODEL
TYPE(FIELDS)       :: YRFIELDS
TYPE(MTRAJ)        :: YRMTRAJ

TYPE(CLASS_VARBC)  :: YVARBC
TYPE(TOVSCV)      :: YLTCV
TYPE(TOVSCV_BGC)  :: YLTCV_BGC
TYPE(JO_TABLE)     :: YLJOT    ! Jo-table
CLASS(DBASE), ALLOCATABLE :: YLODB
LOGICAL :: LLCOUPACTIVE

TYPE(TFPOS) :: YLFPOS
!     ------------------------------------------------------------------

#include "ifs_init.intfb.h"
#include "abor1.intfb.h"
#include "cnt1.intfb.h"
#include "cprep1.intfb.h"
#include "gstats_output_ifs.intfb.h"
#include "su0yoma.intfb.h"
#include "su0yomb.intfb.h"
#include "final_stats.intfb.h"
#include "cseaice.intfb.h"
#include "csstbld.intfb.h"
#include "incli0.intfb.h"
#include "cprep3.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CNT0',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
!*       1.    Initialize level 0 commons.
!              ---------------------------

CALL MODEL_SET(YRMODEL)
#if defined(WITH_OASIS)
IF (PRESENT(LDCOUPACTIVE)) THEN
  LLCOUPACTIVE=LDCOUPACTIVE
ELSE
  LLCOUPACTIVE=.FALSE.
ENDIF
YRMODEL%YRML_AOC%YRMCC%CPLNG_ACTIVE=LLCOUPACTIVE
#endif
!*       1.0   Statistics gathering

IF(LLFIRSTCALL) THEN
  CALL GSTATS(0,0)
  LLFIRSTCALL = .FALSE.
ENDIF

!        1.1   Startup

CALL IFS_INIT(KCOMM=KCOMM)

!     ------------------------------------------------------------------

IF (NCONF /= 903 .AND. NCONF /= 904) THEN

  !*       1.2   Setups up to and including main array allocations

  CALL SU0YOMA(YRGEOMETRY,YRFIELDS%YRSURF,YRMODEL)

  !! point the YGFL pointer in self%yrgfl to that in the model passed in
  YRFIELDS%YRGFL%YGFL => YRMODEL%YRML_GCONF%YGFL
  !! point the model pointer in self to the model passed in
  YRFIELDS%STATE_MODEL => YRMODEL

  !*       1.3   Setups after main array allocations



  CALL SU0YOMB(YLFPOS,YRGEOMETRY,YRFIELDS,YRMTRAJ,YRMODEL,YLJOT,YVARBC,YLTCV,YLTCV_BGC,YLODB)

ELSE

  ! conf 903/904 moves aways from su0yoma/su0yomb
  WRITE(UNIT=NULOUT,FMT='('' POST-PROCESSING JOB'',2I6)') NCONF
  IF (NCONF == 903) THEN
    CALL CPREP3(ITER)
    ! Use the model step to diagnose the computational cost at the end of the job
    NSTEP=ITER
  ELSEIF (NCONF == 904) THEN
    ! test program for change of resolution of an object FIELDS
    !!CALL CPREP4
    CALL ABOR1("CPREP4 is most likely no longer used, and deprecated in preparation for deletion")
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       2.    Decide on configuration of job.
!              -------------------------------

ICONFI=IABS(NCONF)/100

!*       2.0   Integration job.

IF(ICONFI == 0)THEN
  WRITE(UNIT=NULOUT,FMT='('' INTEGRATION JOB'',2I6)')ICONFI,NCONF
  CALL CNT1(YRGEOMETRY,YRFIELDS,YRMTRAJ,YRMODEL,YLJOT,YVARBC,YLTCV,YLODB,YDFPOS=YLFPOS)

!*       2.1   Variational job.

ELSEIF(ICONFI == 1)THEN




  CALL ABOR1('ICONFI==1 without OBS, DA, VARBC')


!*       2.2   2D integration job.

ELSEIF(ICONFI == 2)THEN
  WRITE(UNIT=NULOUT,FMT='('' 2D INTEGR.  JOB'',2I6)')ICONFI,NCONF
  CALL CNT1(YRGEOMETRY,YRFIELDS,YRMTRAJ,YRMODEL,YLJOT,YVARBC,YLTCV,YLODB)

!*       2.3   Kalman Filter / soil moisture analysis

ELSEIF(ICONFI == 3)THEN




    CALL ABOR1('ICONFI==3 without OBS or VARBC')


!*       2.4   Test of the adjoint.

ELSEIF(ICONFI == 4)THEN




    CALL ABOR1('ICONFI==4 without OBS or VARBC')


!*       2.5   Test of the tangent linear.

ELSEIF(ICONFI == 5)THEN




    CALL ABOR1('ICONFI==5 without OBS or VARBC')


!*       2.6   Search for most unstable modes.

ELSEIF(ICONFI == 6)THEN
  WRITE(UNIT=NULOUT,FMT='('' UNSTABLE MODES     '',2I6)')ICONFI,NCONF
  CALL CUN1(YRGEOMETRY,YRFIELDS,YRMTRAJ,YRMODEL,YLJOT,YVARBC)

!*       2.7   Optimal Interpolation job.

ELSEIF(ICONFI == 7)THEN




  CALL ABOR1('ICONFI==7 without ODB')


!*       2.8   Sensitivity job

ELSEIF(ICONFI == 8)THEN
  WRITE(UNIT=NULOUT,FMT='('' SENSITIVITY JOB'',2I6)')ICONFI,NCONF
  CALL CGR1(YRGEOMETRY,YRFIELDS,YRMTRAJ,YRMODEL,YLJOT,YVARBC)

!*       2.9   Setting up the initial conditions.

ELSEIF (NCONF == 923) THEN
  WRITE(UNIT=NULOUT,FMT='('' SETTING UP INITIAL CONDITIONS'',2I6)')ICONFI,NCONF
  CALL INCLI0(YRGEOMETRY,YRFIELDS%YRSURF,YRMODEL%YRML_GCONF%YGFL,&
   & YRMODEL%YRML_PHY_EC%YREPHY,YRMODEL%YRML_PHY_MF,YRMODEL%YRML_AOC%YRMCC)
ELSEIF (NCONF == 931) THEN
  WRITE(UNIT=NULOUT,FMT='('' SETTING UP INITIAL CONDITIONS'',2I6)')ICONFI,NCONF
  CALL CSSTBLD(YRGEOMETRY)
ELSEIF (NCONF == 932) THEN
  CALL CSEAICE(YRGEOMETRY,YRMODEL%YRML_PHY_MF%YRPHY1)
ELSEIF (NCONF == 933) THEN
  IF (MYPROC==1) CALL CSSTBLD(YRGEOMETRY)
  IF(.NOT.LELAM) CALL CSEAICE(YRGEOMETRY,YRMODEL%YRML_PHY_MF%YRPHY1)
ELSEIF (NCONF == 901) THEN
  WRITE(UNIT=NULOUT,FMT='('' SETTING UP INITIAL CONDITIONS'',2I6)')ICONFI,NCONF
  CALL CPREP1(YRGEOMETRY,YRFIELDS%YRSURF,YRMODEL%YRML_PHY_EC%YREPHY,YRMODEL%YRML_GCONF,YRMODEL%YRML_PHY_MF%YRPHY1, &
 &                  YRMODEL%YRML_LBC%TEFRCL)
ELSEIF (NCONF /= 903 .AND. NCONF /= 904) THEN
  WRITE(UNIT=NULOUT,FMT='('' ERROR CNT0  JOB'',2I6)')ICONFI,NCONF
  CALL ABOR1('CALLED IN CNT0')
ENDIF

!     ------------------------------------------------------------------

!*       3.    RUN DOWN
!              --------

!*       3.4   CLOSE FDB AND OTHER POST PROCESSING FILES

IF (LFDBERR .OR. NOUTTYPE==2) THEN
  CALL CLOSE_IOSTREAM(Y_IOSTREAM_FDB)
ENDIF

!*       3.6 STATISTICS

CALL FINAL_STATS('CNT0',NSTEP,YRMODEL%YRML_DYN%YRDYNA%LNHDYN)

IF (LSLDEBUG) CALL UNUSED_HALO_STATS(YRMODEL%YRML_DYN%YRSL, &
 & YRMODEL%YRML_DYN%YRAD, YRMODEL%YRML_PHY_RAD%YRRI, &
 &                                   YRMODEL%YRML_PHY_RAD%YRRO)
CALL IOSTREAM_STATS


CALL GSTATS(0,1)
CALL GSTATS_OUTPUT_IFS(YRMODEL%YRML_GCONF%YRRIP)

WRITE(NULOUT,*) ' *** END CNT0 *** '

IF(.NOT.LMPOFF) THEN

  IF (NCONF /= 903 .AND. NCONF /= 904) THEN
    CALL MTRAJ_DELETE(YRMTRAJ)
    IF (NCONF /= 901) CALL FIELDS_DELETE(YRFIELDS)
    CALL GEOMETRY_DELETE(YRGEOMETRY)
    CALL MODEL_DELETE(YRMODEL)
!!  IF(ALLOCATED(YLODB)) CALL YLODB%DESTROY() -- blows up with Cray CCE 8.6.2 in RAPS17 -- harmless to comment out ?
  ENDIF

!           gather memory usage statistics
  CALL GETMEMSTAT(NULOUT, 'CNT0')

  IF (NPROC > 1) THEN
    IF(LOUTPUT) CALL FLUSH(NULOUT)
    CALL MPL_BARRIER(CDSTRING='CNT0:')
  ENDIF
ENDIF

#ifdef WITH_ATLAS
CALL ATLAS_LIBRARY%FINALISE()
#endif
!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CNT0',1,ZHOOK_HANDLE)
END SUBROUTINE CNT0
