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

SUBROUTINE IFS_INIT(CDEXPV,KCOMM)

! Initialization of IFS, invoked from both OOPS and IFS control flow.
! Initialises global modules, constants etc. and message passing


USE PARKIND1,   ONLY : JPRB, JPIM
USE YOMHOOK,    ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0,     ONLY : LELAM, LGRIB_API, LOBS
USE YOMCT3,     ONLY : NSTEP
USE YOMLUN,     ONLY : NULNAM, NULOUT
USE YOMMP0,     ONLY : MYPROC, NPROC, NPRINTLEV,LOUTPUT
USE YOMNUD,     ONLY : LNUDG
USE YOMRIP0,    ONLY : NINDAT, NSSSSS
USE YOMVAR,     ONLY : SETUP_VAR
USE YOMCST,     ONLY : SETUP_CONSTANTS
USE YOMOBS,     ONLY : LVDFTRAJ
USE GRIB_HANDLES_MOD, ONLY : READ_GRIB_SAMPLES
USE YOMARG,     ONLY : CNMEXP
USE INTDYN_MOD, ONLY : SUINTDYN
USE YOMDYNA_STATIC,ONLY : SUDYNA_STATIC
USE ENKF_MIX  , ONLY : SETUP_ENKF
USE SATS_MIX  , ONLY : SUSATS
USE MPL_MODULE, ONLY : MPL_NPROC, MPL_MYRANK
USE STACK_MIX  ,ONLY : INIT_STACK
USE MPL_MPIF   ,ONLY : MPI_COMM_WORLD
USE MPL_DATA_MODULE, ONLY : MPL_COMM
USE FIELD_DEFINITIONS, ONLY : READ_DYNAMIC_NAMESPACE
#ifdef WITH_FCKIT
USE FCKIT_MODULE, ONLY : FCKIT_MAIN, FCKIT_LOG, FCKIT_EXCEPTION, FCKIT_EXCEPTION_HANDLER
#endif
#ifdef WITH_ATLAS
USE ATLAS_MODULE       , ONLY : ATLAS_LIBRARY
#endif

IMPLICIT NONE

CHARACTER (LEN=4),  INTENT(IN), OPTIONAL :: CDEXPV
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KCOMM

CHARACTER (LEN=35) ::  CLINE
CHARACTER(LEN=1)   :: CLEC_MEMINFO
REAL(KIND=JPHOOK ) :: ZHOOK_HANDLE
#ifdef WITH_FCKIT
PROCEDURE(FCKIT_EXCEPTION_HANDLER), POINTER :: FUNPTR
#include "include_abor1_intfb.h"
#endif

#include "ec_meminfo.intfb.h"
#include "suarg.intfb.h"
#include "suct0.intfb.h"
#include "suect0.intfb.h"
#include "suct1.intfb.h"
#include "sucst_ifsaux.h"
#include "sudyncore.intfb.h"
#include "sujfh.intfb.h"
#include "sulun.intfb.h"
#include "sump0.intfb.h"
#include "sumpini.intfb.h"
#include "sumpini_prt.intfb.h"
#include "sumsc.intfb.h"
#include "suoph0.intfb.h"
#include "surip0.intfb.h"
#include "sutrans0.intfb.h"
#include "suetrans0.intfb.h"
#include "sufp_ctl.intfb.h"
#include "suppvi.intfb.h"
#include "sutim.intfb.h"
#include "suini.intfb.h"
IF (LHOOK) CALL DR_HOOK('IFS_INIT',0,ZHOOK_HANDLE)

CLINE='----------------------------------'

IF (PRESENT(CDEXPV)) CNMEXP = CDEXPV

!*    Determine if message passing or shared memory version is used
CALL SUMPINI

NPROC=MPL_NPROC()
MYPROC=MPL_MYRANK()
! Every program needs to be initialised
#ifdef WITH_FCKIT
CALL FCKIT_MAIN%INITIALISE()

! Register ABOR1 as fckit's exception handler
FUNPTR => ABOR1_EXCEPTION_HANDLER
CALL FCKIT_EXCEPTION%SET_HANDLER( FUNPTR )
#endif

CALL INIT_STACK(1)
IF (PRESENT(KCOMM)) THEN
  CALL EC_MEMINFO(-1,"oops:ifs_init",KCOMM,KBARR=1,KIOTASK=0,KCALL=0)
ELSE
  CALL GET_ENVIRONMENT_VARIABLE('EC_MEMINFO',CLEC_MEMINFO)
  IF (CLEC_MEMINFO /= '0') CALL ABOR1("IFS_INIT needs an MPI communicator for call to EC_MEMINFO")
ENDIF

!*    Initialize YOMLUN (sets up NULOUT and other constatns)
CALL SULUN
#ifdef WITH_FCKIT
IF(LOUTPUT) THEN
  CALL FCKIT_LOG%SET_FORTRAN_UNIT(UNIT=NULOUT,STYLE=FCKIT_LOG%PREFIX)
  IF(MYPROC == 1) CALL FCKIT_LOG%ADD_STDOUT(STYLE=FCKIT_LOG%PREFIX)
ELSE
  CALL FCKIT_LOG%RESET()
ENDIF
CALL FCKIT_LOG%INFO('LOGGER RE-CONFIGURED BY IFS_INIT')
#endif
WRITE(NULOUT,*) '--- End of logical units setup ------',CLINE

!*    Prints quantities computed in SUMPINI (printings must be done after calling SULUN).
CALL SUMPINI_PRT

!*    Initialize machine-specific constants
WRITE(NULOUT,*) '--- Set up machine-specific constants',CLINE
CALL SUMSC(NULOUT)

!*   (Start) Initialize Universal constants
CALL SUCST_IFSAUX

!*    Get command line arguments
WRITE(NULOUT,*) '--- Set up command line arguments ---',CLINE
CALL SUARG

!*    Initialize level 0 control common
WRITE(NULOUT,*) '--- Set up control common 0 ---------',CLINE
CALL SUCT0(NULOUT)
IF (LELAM)  CALL SUECT0(NULOUT)

!*    Initialize level 1 control common
WRITE(NULOUT,*) '--- Set up control common 1 ---------',CLINE
CALL SUCT1

!*    Initialize full post processing module (data shared by all objects)
WRITE(NULOUT,*) '------ Set up F-post processing module ',CLINE
CALL SUFP_CTL

!*    Initialize vertical interpolator (used both by obs operators and post-processing)
WRITE(NULOUT,*) '---- Set up vertical interpolator ---',CLINE
CALL SUPPVI

! YOMCT3
NSTEP      = 0

! YOMOBS
LVDFTRAJ = .FALSE.

! Initialize MASS VF Option
WRITE(NULOUT,*) '--- Set up MASS VF Option ---------',CLINE
CALL SUJFH

!*    Initialize dynamics: part A
WRITE(NULOUT,*) '--- Set up dynamics part A static ---------',CLINE
CALL SUDYNA_STATIC

!*      YOMCVER setup (VFE keys)
! Moved to GEOMETRY

!*    Initialize message passing interface
WRITE(NULOUT,*) '--- Set up message passing interface ',CLINE
CALL SUMP0(NULOUT,NULNAM)

!*    Initialize control of variational assimilation
WRITE(NULOUT,*) '---- Set up variational assimilation ',CLINE
CALL SETUP_VAR(NULOUT)
!*    DISABLED: Initialize nudging
LNUDG = .FALSE.
!WRITE(NULOUT,*) '------ Set up nudging ',CLINE
!CALL SUNUD(NULOUT)

!*       dynamical core setup
WRITE(NULOUT,*) '------ Set up dynamical core -------------',CLINE
CALL SUDYNCORE

!*    Initialize constants
WRITE(NULOUT,*) '------ Set up constants -------------',CLINE
CALL SETUP_CONSTANTS(NULOUT,NPRINTLEV)

!*    Initialize model time
WRITE(NULOUT,*) '------ Set up model time (YOMRIP0) --',CLINE
CALL SURIP0(NULOUT)  ! all but LASTRF will be over-written

!*    Initialize control of DFI
WRITE(NULOUT,*) '------ Set up DFI initialization ',CLINE
CALL SUINI

!*    Initialize file handling
WRITE(NULOUT,*) '---- Set up files : names, FA --',CLINE
IF (PRESENT(CDEXPV)) THEN
  CALL SUOPH0(CDEXPV)
!ELSE
! already called by suarg with cnmexp
!  CALL SUOPH0(CNMEXP)
ENDIF

!*    Initialize transform package
WRITE(NULOUT,*) '------ Initialize transform package --',CLINE
IF (LELAM) THEN
  CALL SUETRANS0('SETUP')
ELSE
  CALL SUTRANS0('SETUP')
ENDIF

!* Initialize YOMRINC (bit odd here but ...) MH




! Read GRIB API samples (templates)

IF(LGRIB_API) THEN
  CALL READ_GRIB_SAMPLES
ENDIF

! Read dynamic field names for field_containers
WRITE(NULOUT,*) '------ Read dynamic namespace ---------',CLINE
CALL READ_DYNAMIC_NAMESPACE

!*    Initialize some structures used in the adiabatic model or some GP.. routines
WRITE(NULOUT,*) '---- Set up some structures used in the adiabatic model -------'
CALL SUINTDYN

!*    Initialize computer time
WRITE(NULOUT,*) '------ Set up computer time ---------',CLINE
CALL SUTIM(NULOUT)

!*    Initialize control of different TL/AD diagnostics




!*     Initialize Atlas
#ifdef WITH_ATLAS
CALL ATLAS_LIBRARY%INITIALISE( COMM=MPL_COMM )
#endif

! MC - Do we need it?  MH Yes!
!*    Initialize Ensemble Kalman Filter options
WRITE(NULOUT,*) '--- Set up ENKF options ---------',CLINE
CALL SETUP_ENKF


!*    Initialize Satellites sensors
IF (.NOT.LOBS) THEN
  WRITE(NULOUT,*) '------ Set up satellites sensors -------',CLINE
  CALL SUSATS
! For now, if (LOBS) SUSATS is still called in GETSATID
ENDIF

IF (LHOOK) CALL DR_HOOK('IFS_INIT',1,ZHOOK_HANDLE)

END SUBROUTINE IFS_INIT
