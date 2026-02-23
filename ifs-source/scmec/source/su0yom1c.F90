! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU0YOM1C(YDGEOMETRY,YDMODEL,YDFIELDS)

USE PARKIND1,           ONLY : JPIM, JPRB
USE GEOMETRY_MOD,       ONLY : GEOMETRY
USE TYPE_MODEL,         ONLY : MODEL
USE FIELDS_MOD,         ONLY : FIELDS
USE YOMHOOK,            ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN,             ONLY : NULOUT, NULNAM
USE YOMRIP0,            ONLY : NINDAT, NSSSSS
USE YOMMP0,             ONLY : LSCMEC, MYPROC, LOUTPUT, NPRINTLEV
USE YOMLOG1C,           ONLY : OUTFORM
USE YOMCST,             ONLY : SETUP_CONSTANTS
USE INTDYN_MOD,         ONLY : SUINTDYN
USE SURFACE_FIELDS_MIX, ONLY : ALLO_SURF
USE VARIABLES_MOD,      ONLY : VARIABLES, VARIABLES_CREATE, VARIABLES_DELETE

#ifdef DOC

!**** *SU0YOM1C*  - INITIALIZE LEVEL 0 COMMONS

!     PURPOSE.
!     --------
!           INITIALIZE LEVEL 0 COMMONS (CONSTANT ALONG ALL THE JOB).

!**   INTERFACE.
!     ----------
!        *CALL* *SU0YOM1C*

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
!      SULUN     - INITIALIZE LOGICAL UNITS
!      SUARG1C   - GET COMMAND LINE ARGUMENTS
!      SUCT0     - INITIALIZE CONTROL COMMON BLOCK OF LEVEL 0
!      SURIP1C   - INITIALIZE MODEL TIME
!      SUCST     - INITIALIZE CONSTANTS
!      SU0PHY    - INITIALIZE PHYSICS SWITCHES
!      SUGFL1C   - INITIALIZE GFL FIELDS
!      SUNDDH    - INITIALIZE DDH
!      SUALLO1C  - ALLOCATE SPACE FOR ARRAYS
!      SUDYN1C   - INITIALIZE DYNAMICS
!      SUPHEC    - INITIALIZE PHYSICS
!      SU1C      - INITIALIZE ONE-COLUMN SWITCHES

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the one-column model

!     AUTHOR.
!     -------
!        Joao Teixeira and Pedro Viterbo  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL      93-12-15
!        Joao Teixeira  Nov.96   all routines called are ...1C +
!                                calling suallo1c (f90).
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!     ------------------------------------------------------------------
#endif


IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
TYPE(FIELDS),   INTENT(INOUT) :: YDFIELDS

INTEGER(KIND=JPIM) :: IGFLCONF
TYPE(VARIABLES)    :: YNLVARS

REAL(KIND=JPHOOk)  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "sup1c_nc.intfb.h"
#include "sudyn1c.intfb.h"
#include "sumcc.intfb.h"
#include "sudyncore.intfb.h"
#include "sudyna.intfb.h"
#include "sulun.intfb.h"
#include "sunddh.intfb.h"
#include "suct0.intfb.h"
#include "sujfh.intfb.h"
#include "su0phy.intfb.h"
#include "suallo1c.intfb.h"
#include "suphec.intfb.h"
#include "suphyds1c.intfb.h"
#include "sufa.intfb.h"
#include "sud1c_nc.intfb.h"
#include "su1c.intfb.h"
#include "sugmv1c.intfb.h"
#include "sugfl1c.intfb.h"
#include "surip.intfb.h"
#include "surip1c.intfb.h"
#include "suarg1c.intfb.h"
#include "su_surf_flds.intfb.h"
#include "susavtend.intfb.h"
#include "sugeometry1c.intfb.h"
#include "sudimf1.intfb.h"
#include "sudimf2.intfb.h"
#include "sugfl.intfb.h"
#include "susc2c.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0YOM1C',0,ZHOOK_HANDLE)

!------------------------------------------------------------------
!        CALLS FROM su0yoma.F90
!------------------------------------------------------------------

!        1.    INITIALIZE LOGICAL UNITS AND COMMAND LINE ARGUMENTS.
!              ----------------------------------------------------

LSCMEC = .TRUE.
!hack, I know...
YDGEOMETRY%YRDIM%NGPBLKS=1
YDGEOMETRY%YRGEM%NGPTOT=1
YDGEOMETRY%YRGEM%NGPTOTG=1

!        and also MYPROC=1
MYPROC=1
LOUTPUT = .TRUE. ! Set logical to turn PE0 output on
NPRINTLEV = 0    ! Set level of print output

CALL SULUN
CALL SUARG1C

!        2.    INITIALIZE LEVEL 0 CONTROL COMMON.
!              ----------------------------------

CALL SUCT0(NULOUT)

!*    Setup model geometry
CALL SUGEOMETRY1C(YDGEOMETRY,YDMODEL%YRML_DYN)

!     Initialize MASS VF Option
CALL SUJFH

!        3.    INITIALIZE DYNAMICAL CORE (from sudyna.F90) .
!              ---------------------------------------------

CALL SUDYNCORE
CALL SUDYNA(YDGEOMETRY%YRDIM,YDMODEL%YRML_DYN%YRDYNA,YDGEOMETRY%YRCVER%LVERTFE, &
 & YDGEOMETRY%YRCVER%NDLNPR,YDGEOMETRY%LNONHYD_GEOM,NULOUT)


!        4.    INITIALIZE MODEL TIME.
!              ----------------------

CALL SURIP(YDGEOMETRY%YRDIM,YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_GCONF%YRRIP)
CALL SURIP1C(NULOUT)

!        5.    INITIALIZE CONSTANTS.
!              ---------------------

!CALL SETUP_CONSTANTS(NULOUT,NINDAT)

!        6.    INITIALIZE PHYSICS SWITCHES.
!              ----------------------------

CALL SU0PHY(YDMODEL,NULOUT)

!        7.    INITIALIZE PHYSICS FIELDS DESCRIPTORS.
!              --------------------------------------

CALL SUPHYDS1C(YDMODEL)
CALL SUFA

!*       8. GEOMETRY SET-UP.
!           ----------------

!*    Initialize some structures used in the adiabatic model or some GP.. routines
!     Must be called before calling SUSTA and any routine calling GPHPRE.
CALL SUINTDYN


!*    Initialize number of fields dimensions (part 1)
CALL SUDIMF1(YDMODEL)
CALL SUDIMF2(YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYNA)

!*    Initialize pointers to save ECMWF physical tendencies

CALL SUSAVTEND(YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_EC%YREPHY%LSLPHY,YDMODEL%YRML_PHY_G%YRSLPHY)

!        9.    SET UP UNIFIED TREATMENT OF GRID-POINT FIELDS.
!              ----------------------------------------------

CALL SUGMV1C(YDFIELDS%YRGMV)

IGFLCONF=1
CALL SUGFL(YDGEOMETRY%YRDIMV,YDMODEL,IGFLCONF)

!       10.    SET UP SURFACE GRID-POINT FIELDS.
!              ---------------------------------

CALL SU_SURF_FLDS(YDGEOMETRY%YRDIMV,YDFIELDS%YRSURF,YDMODEL)


!       11.    INITIALIZE DDH (horizontal domains diagnostics).
!              ------------------------------------------------
!              (partially, for callpar consistent call)

CALL SUNDDH(YDGEOMETRY,YDMODEL)

!       12.    ALLOCATION OF POINTERS.
!              -----------------------

CALL SUALLO1C(YDGEOMETRY,YDMODEL)

!       12.1   SURFACE FIELDS

!CALL ALLO_SURF(YDGEOMETRY%YRDIM,YDFIELDS%YRSURF)



!------------------------------------------------------------------
!        CALLS FROM su0yomb.F90
!------------------------------------------------------------------

!       13.    INITIALIZE GEOMETRY.
!              --------------------
!*    Initialize some structures used in the adiabatic model or some GP.. routines
!     Must be called before calling SUSTA and any routine calling GPHPRE.
CALL SUINTDYN


! Instead of suorog, allocate the space for grid-point orography fields
YDGEOMETRY%YRDIM%NGPBLKS=1
ALLOCATE(YDGEOMETRY%YROROG(YDGEOMETRY%YRDIM%NGPBLKS))
ALLOCATE(YDGEOMETRY%YROROG(YDGEOMETRY%YRDIM%NGPBLKS)%OROG (YDGEOMETRY%YRDIM%NPROMA))
ALLOCATE(YDGEOMETRY%YROROG(YDGEOMETRY%YRDIM%NGPBLKS)%OROGL(YDGEOMETRY%YRDIM%NPROMA))
ALLOCATE(YDGEOMETRY%YROROG(YDGEOMETRY%YRDIM%NGPBLKS)%OROGM(YDGEOMETRY%YRDIM%NPROMA))

CALL VARIABLES_CREATE(YNLVARS, .FALSE.)
CALL SUSC2C(YDGEOMETRY,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_GCONF,YDMODEL%YRML_PHY_MF%YRPHY, &
 &          YDMODEL%YRML_DYN%YRDYNA,YDMODEL%YRML_LBC%LTENC,  YNLVARS,YDFIELDS%YRGFL, &
 &          YDFIELDS%YRGMV,YDFIELDS%YRSURF)
CALL VARIABLES_DELETE(YNLVARS)

!       14.    INITIALIZE DYNAMICS.
!              -------------------

CALL SUDYN1C(YDGEOMETRY,YDMODEL)

! ---- Set up MCC climate model keys ----
CALL SUMCC(YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_EC%YREPHY,NULOUT)

!       15.    INITIALIZE PHYSICS.
!              -------------------

CALL SUPHEC(YDGEOMETRY,YDMODEL,NULOUT)

!       16.    INITIALIZE STOCHASTIC PHYSICS.
!              ------------------------------

! if this is needed only copy the relevant part here from STOPH_MIX explicitly
!CALL SURAND1

!------------------------------------------------------------------
!        CALLS for scm specifics
!------------------------------------------------------------------

!       17.    INITIALIZE ONE-COLUMN SWITCHES.
!              -------------------------------

CALL SU1C(YDMODEL)

!       18.    Open and set up NetCDF datafiles.
!              ---------------------------------

IF (OUTFORM .EQ. 'netcdf') THEN
  CALL SUP1C_NC(YDGEOMETRY%YRDIMV,YDMODEL,YDFIELDS%YRSURF)
  CALL SUD1C_NC(YDGEOMETRY%YRDIMV,YDMODEL,YDFIELDS%YRSURF)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0YOM1C',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SU0YOM1C
