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

SUBROUTINE SU0YOMA(YDGEOMETRY,YDSURF,YDMODEL)

!**** *SU0YOMA*  - INITIALIZE LEVEL 0 COMMONS AND SOME HIGHER (PART 1)

!     PURPOSE.
!     --------
!           INITIALIZE LEVEL 0 COMMONS (CONSTANT ALONG ALL THE JOB)
!        AND SOME HIGHER LEVEL COMMONS. CALLS ROUTINES THAT PERFORM
!        ALL THE PREPARATIONS NEEDED TO EXECUTE MODEL. THE TASK OF
!        INITIALIZING THE COMMONS IS DIVIDED BETWEEN TWO ROUTINES
!        (SU0YOMA AND SU0YOMB) TO AVOID PROBLEMS WITH USING POINTER
!        ARRAYS WHOSE DIMENSIONS ARE NOT DEFINED UNTIL AFTER CALLING SUDIM.

!**   INTERFACE.
!     ----------
!        *CALL* *SU0YOMA*

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
!      MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!      ORIGINAL      :87-10-15

!     MODIFICATIONS.
!     --------------
!      R. El Khatib  : 01-01-26 SUCDDH moved inside SUNDDH
!      R. El Khatib  : 01-02-01 SUCFU/XFU moved to su0yomb after call susta
!      Modified 08-2002 C. Smith : use "w" as prognostic variable in the
!       semi-lag advection of vertical divergence in the NH model.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      C. Fischer    : 04-02-25 Merge Aladin into Arpege/IFS su0yoma
!      T.Kovacic     :04-03-01 Diagnostics on physical tendencies setup
!      R. El Khatib  : 04-08-05 suedim moved inside sudim2
!      15-Mar-2005 R. El Khatib SUJFH
!      C. Temperton  : 06-09-05 avoid calling SUJB if NCONF=201
!      A. Fouilloux  : do not set up arrays for observations
!      R. El Khatib  24-Jul-2009 Setup satellite sensors if .not.LOBS
!      K. Yessad (Aug 2009): prune conf 912, externalise conf 911.
!      R. El Khatib  29-Mar-2010 suclmicst moved in the setup.
!      K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!      M. Fisher   7-March-2012 Move Jb setup into SU0YOMB
!      K. Yessad (Nov 2011): call to SUSAVTEND moved from SUDYN to SU0YOMA.
!      B. Bochenek (Apr 2012): call to SUNDDH moved to SU0YOMB
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      W. Deconinck (Feb 2015): Introduce Atlas framework
!      R. El Khatib : 01-Dec-2014 move up sutim
!      K. Yessad (Oct 2015): move SU0PHY+SUTRAJP+SUNUD just after SUGEOMETRY.
!      M. Janiskova 11-Apr-2016: Call TL/AD diagnostic setup
!      R. El Khatib 17-Aug-2016 move suoph up from su0yomb to su0yoma and move down sufa
!      S. Massart   19-Feb-2019 Augmented control variable
!      F. Vana      11-Sep-2020 Cleaning & moving SLAVEPP setup to better location
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMCT0             , ONLY : LMINIM, LR2D, NCONF, LOBSC1
USE YOMCLMICST         , ONLY : SETUP_CLMICST

USE YOMSATSIM,    ONLY : NSATSIM


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TSURF)   , INTENT(INOUT) :: YDSURF
TYPE(MODEL)   , INTENT(INOUT) :: YDMODEL
CHARACTER (LEN = 35) ::  CLINE
INTEGER(KIND=JPIM) :: IGFLCONF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "su0phy.intfb.h"
#include "suallo.intfb.h"
#include "suchet.intfb.h"
#include "sudimf1.intfb.h"
#include "sudimf2.intfb.h"
#include "sudim_traj.intfb.h"
#include "sudyna.intfb.h"




#include "suechk.intfb.h"
#include "sufa.intfb.h"
#include "sugeometry.intfb.h"
#include "sugfl.intfb.h"
#include "sumcc.intfb.h"
#include "sunud.intfb.h"
#include "surip.intfb.h"



#include "sutrajp.intfb.h"
#include "suoph.intfb.h"
#include "su_surf_flds.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0YOMA',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

CLINE='----------------------------------'

CALL GSTATS(1934,0)

!     ------------------------------------------------------------------

!*       3. GEOMETRY SET-UP.
!           ----------------

!*    Setup model geometry
WRITE(NULOUT,*) '------ Set up model geometry  ------',CLINE
CALL SUGEOMETRY(YDGEOMETRY)

!     ------------------------------------------------------------------
!*       4. "MODEL" PART SET-UP.
!           --------------------
WRITE(NULOUT,*) '--- Set up dynamics part A ---------',CLINE
CALL SUDYNA(YDGEOMETRY%YRDIM,YDMODEL%YRML_DYN%YRDYNA,YDGEOMETRY%YRCVER%LVERTFE, &
 & YDGEOMETRY%YRCVER%NDLNPR,YDGEOMETRY%LNONHYD_GEOM,NULOUT)

!*    Initialize YOMRIP variables
WRITE(NULOUT,*) '------ Set up YOMRIP variables ',CLINE
CALL SURIP(YDGEOMETRY%YRDIM,YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_GCONF%YRRIP)

!*    Initialize control of physical parameterizations
WRITE(NULOUT,*) '-- Set up physical parameterizations ',CLINE
CALL SU0PHY(YDMODEL,NULOUT)

!*    Initialize control of the adjoint physics
WRITE(NULOUT,*) '------ Set up NH trajectory -------',CLINE
CALL SUTRAJP

!*    Initialize nudging
WRITE(NULOUT,*) '------ Set up nudging ',CLINE
CALL SUNUD(YDMODEL%YRML_PHY_MF%YRARPHY,NULOUT)

!*    Initialize special keys for the climate version
IF (.NOT.LR2D) THEN
  WRITE(NULOUT,*) '---- Set up MCC climate model keys --',CLINE
  CALL SUMCC(YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_EC%YREPHY,NULOUT)
ENDIF

!*    Initialize some dimensions for trajectory and background.
WRITE(NULOUT,*) '------ Set up some dimensions for trajectory and background ------',CLINE
CALL SUDIM_TRAJ(YDGEOMETRY%YRDIM)

!*    Initialize file handling
WRITE(NULOUT,*) '---- Set up files handling, FA --',CLINE
CALL SUOPH(YDGEOMETRY)

!*    Initialize Arpege field names and GRIB packing options
WRITE(NULOUT,*) '--- Set up GRIB packing options ---',CLINE
CALL SUFA

!*    Initialize number of fields dimensions (part 1)
WRITE(NULOUT,*) '------ Set up number of fields dimensions, part 1  ------',CLINE
CALL SUDIMF1(YDMODEL)

!*    Set up unified_treatment grid-point fields
WRITE(NULOUT,*) '---- Set up unified_treatment grid-point fields -',CLINE
IF (ANY(SPREAD(NCONF,1,5) == (/1,302,601,801,903/))) THEN
  IGFLCONF=1
ELSE
  IGFLCONF=NCONF
ENDIF
CALL SUGFL(YDGEOMETRY%YRDIMV,YDMODEL,IGFLCONF)

!*    Enable simulated satellite images in forecast mode only
IF (NCONF == 1 .AND. .NOT. LOBSC1) THEN




  NSATSIM=0

ENDIF

!*    Set up for surface grid-point fields
WRITE(NULOUT,*) '---- Set up for surface grid-point fields ----',CLINE
CALL SU_SURF_FLDS(YDGEOMETRY%YRDIMV,YDSURF,YDMODEL)

!*    Initialize number of fields dimensions (part 2)
WRITE(NULOUT,*) '------ Set up number of fields dimensions, part 2  ------',CLINE
CALL SUDIMF2(YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA)

!*    Initialize minimisation I/O scheme if necessary
IF (LMINIM) THEN
  WRITE(NULOUT,*) '------ Set up minimization I/O -------',CLINE
  CALL SUIOMI(NULOUT)
ENDIF

!*    Initialize gridpoint evolution diagnostics
WRITE(NULOUT,*) '------ Set up gridpoint diagnostics ----',CLINE
CALL SUECHK(YDGEOMETRY,YDMODEL%YRML_GCONF%YRDIMF)

!*    Initialize extended control variable options
WRITE(NULOUT,*) '--- Set up ECV geometry ---------',CLINE




!*    Initialize diagnostics on physical tendencies
WRITE(NULOUT,*) '------ Set up diagnostics on phys. tendencies -',CLINE
CALL SUCHET(YDGEOMETRY%YRDIMV)

!*    Initialise constants of microphysics scheme ICE3
WRITE(NULOUT,*) '------ Set up  constants of microphysics scheme ICE3 -',CLINE
CALL SETUP_CLMICST

!*    Allocate grid point and spectral arrays
WRITE(NULOUT,*) '------ Set up : array allocations ------',CLINE
CALL SUALLO(YDGEOMETRY,YDMODEL)

CALL GSTATS(1934,1)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0YOMA',1,ZHOOK_HANDLE)
END SUBROUTINE SU0YOMA
