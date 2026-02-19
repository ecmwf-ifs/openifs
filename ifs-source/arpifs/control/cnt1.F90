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

#define __FILENAME__ "cnt1.F90"
SUBROUTINE CNT1(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDODB,YDFPOS)

!**** *CNT1*  - Controls integration job at level 1.

!     Purpose.
!     --------
!           Controls the integration job at level 1.

!**   Interface.
!     ----------
!        *CALL* *CNT1

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
!        SU1YOM  - initialize level 1 commons
!        SUOBS   - initialize obs processing
!        CNT2    - control integration job at level 2
!        SCREEN  - observation screening
!        SACMAC1 - save CMA files and evaluate cost function
!        PERTOBS - perturb observations

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      R. El Khatib: 02-11-12  Pruning LREFFP
!      D. Dee     : 04-03-03  Variational bias correction
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      M.Fisher : 21-10-04  Perturb observations
!      P.Burton :  29-09-04 P.Burton LVARBC cleanup
!      Dick Dee : 01-03-08  New VarBC setup
!      Y.Tremolet    17-Nov-2008 Jc-DFI diagnostics in outer loop
!      R. El Khatib: 09-08-28  Disable call to  deallocate subroutines  when NFPCT0>0
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      A. Geer      27-Jul-2015   VarBC is now an object, for OOPS
!      A. Geer      11-Feb-2016 Supergom object for model<->observation space interpolation
!      P. Lean      16-Aug-2016 ODB is now an object, for OOPS
!      P. Lean      22-Mar-2017 Jo-table is now an object, for OOPS
!      S. Massart   19-Feb-2019 Paramater optimisation
!     ------------------------------------------------------------------

USE TYPE_MODEL    , ONLY : MODEL
USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE FIELDS_MOD    , ONLY : FIELDS
USE MTRAJ_MOD     , ONLY : MTRAJ
USE PARKIND1      , ONLY : JPRB
USE YOMHOOK       , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE JO_TABLE_MOD  , ONLY : JO_TABLE
USE YOMCT0        , ONLY : LOBS, LOBSC1, LSCREEN, L_SCREEN_CALL
USE YOMVAR        , ONLY : LVARBC, LJCDFI, LECV
USE YOMDIMO       , ONLY : NOBTOT, NACTIM
USE YOMSCC        , ONLY : LPERTURB
USE YOMLUN        , ONLY : NULOUT
USE VARBC_CLASS   , ONLY : CLASS_VARBC
USE TOVSCV_MOD    , ONLY : TOVSCV
USE SUPERGOM_CLASS, ONLY : CLASS_SUPERGOM
USE YOMVRTL       , ONLY : LOBSTL



USE YOMLOCS       , ONLY : TOBSLOCS

USE DBASE_MOD     , ONLY : DBASE
USE FULLPOS       , ONLY : TFPOS
USE YOMJBECV      , ONLY : LECPHYSPARECV, READ_FG_ECV, SAVE_FG_ECV, YRECV5
USE YOMJBPAR1DECV , ONLY : SUPARECVTRAJ, PARECV_SAVE
USE YOMJBECPHYSECV , ONLY : LSOLARCST
#ifdef WITH_ATLAS
USE ATLAS_MODULE  , ONLY : ATLAS_TRACE
#endif

IMPLICIT NONE

TYPE(GEOMETRY)   ,INTENT(INOUT) :: YDGEOMETRY !! INOUT needed for call to CNT2
TYPE(FIELDS)     ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)      ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)      ,INTENT(INOUT) :: YDMODEL
TYPE(JO_TABLE)   ,INTENT(INOUT) :: YDJOT
TYPE(CLASS_VARBC),INTENT(INOUT) :: YDVARBC
TYPE(TOVSCV)     ,INTENT(IN)    :: YDTCV
CLASS(DBASE)     ,INTENT(INOUT) :: YDODB
TYPE(TFPOS)      ,INTENT(IN), OPTIONAL :: YDFPOS

TYPE(CLASS_SUPERGOM) :: YGOM
TYPE(CLASS_SUPERGOM) :: YGOM5
TYPE(TOBSLOCS) :: YOBSLOCS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef WITH_ATLAS
TYPE(ATLAS_TRACE) :: TRACE
#endif

#include "cnt2.intfb.h"
#include "dealxmo.intfb.h"
!      -----------------------------------------------------------------

!*       1.    Initialize YOMCT1.
!              ------------------

IF (LHOOK) CALL DR_HOOK('CNT1',0,ZHOOK_HANDLE)
#ifdef WITH_ATLAS
TRACE = ATLAS_TRACE("cnt1.F90",__LINE__,"CNT1")
#endif





!      -----------------------------------------------------------------
!*       2.    set up observations (LOBSC1)
!              ----------------------------

IF (LOBSC1 .AND. LOBS) THEN
  call abor1("OIFS - cnt1 call to obs setup (LOBSC1 .AND. LOBS) should never be called - EXIT")


ENDIF

IF (LJCDFI) THEN




  call abor1("OIFS - cnt1 call to obs setup (LJCDFI) should never be called - EXIT")

ENDIF

!      -----------------------------------------------------------------
!*       2.5   set up VarBC (LOBSC1)
!              ---------------------

IF(LOBSC1 .AND. LVARBC) THEN



  call abor1("OIFS - cnt1 call to varbc setup (LOBSC1 .AND. LVARBC) should never be called - EXIT")

ENDIF

!      -----------------------------------------------------------------
!*       2.6   set up parameters optimization
!              --------------------------------------
!      -----------------------------------------------------------------

!*       3.    Call level 2 control routine.
!              -----------------------------

CALL CNT2(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YGOM5,YDODB,YDFPOS=YDFPOS)
CALL DEALXMO(YDFIELDS%YRGFL,YDFIELDS%YRGMV)

! Free as much memory as possible
IF (LSCREEN) THEN



  call abor1("OIFS - cnt1 call to obs setup (LSCREEN) should never be called - EXIT")

ENDIF

!      -----------------------------------------------------------------
!*       3.5   perturb observations (LOBSC1)
!              ----------------------------

IF (LOBSC1 .AND. LPERTURB) THEN



  call abor1("OIFS - cnt1 call to varbc setup (LOBSC1 .AND. LPERTURB) should never be called - EXIT")

ENDIF

!      -----------------------------------------------------------------

!*       4.    Observation screening (LSCREEN)
!              -------------------------------

IF (LSCREEN.AND.L_SCREEN_CALL) THEN





  call abor1("OIFS - cnt1 call to obs setup (LSCREEN.AND.L_SCREEN_CALL) should never be called - EXIT")

ENDIF

!      -----------------------------------------------------------------

!*       5.    Save VarBC information for cycling
!              ----------------------------------

IF (LOBS .AND. LVARBC) THEN



  call abor1("OIFS - cnt1 call to varbc setup (LOBS .AND. LVARBC) should never be called - EXIT")

ENDIF

!      -----------------------------------------------------------------

!*       6.    save CMA files and evaluate cost function (LOBSC1)
!              --------------------------------------------------
!      -----------------------------------------------------------------
!*       7.    Save parameters optim. information for cycling
!              ------------------------------------------------------

IF (LOBS .AND. LECV) THEN




  call abor1("OIFS - cnt1 call to obs setup (LSCREEN.AND.L_SCREEN_CALL) should never be called - EXIT")
!  IF (LECPHYSPARECV) CALL PARECV_SAVE()

ENDIF

!      -----------------------------------------------------------------

#ifdef WITH_ATLAS
CALL TRACE%FINAL()
#endif

IF (LHOOK) CALL DR_HOOK('CNT1',1,ZHOOK_HANDLE)
END SUBROUTINE CNT1
#undef __FILENAME__
