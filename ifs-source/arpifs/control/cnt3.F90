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

SUBROUTINE CNT3(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDGOM5,YDODB,YDFPOS)

!**** *CNT3*  - Controls integration job at level 3

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *CNT3

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
!      See includes below.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      Modified : 02-09-30 V.Guidard&C.Fischer - 3dfgat coupling switched on
!      J. Masek : 12-10-2002: Call of SUPONG moved from SUDYN.
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 02-11-12  Pruning LREFFP
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      C. Fischer 04-02-26 Merge Aladin into Arpege/IFS cnt3
!      C. Fischer 04-10-20 call sueqlimsat
!      Y.Tremolet    21-Jul-2004 Model error
!      M. Jidane : 13-04-2006 : SWAP37 NO MORE IN USE
!      M. Drusch:    17-Jan-2007 introduce nconf 302
!      B. Chapnik:  22-sep-2008 allows upspec for aladin
!      G. Desroziers 22-Dec-2008: Enable transf. of ARPEGE file in GRIB format (to be used in femars)
!      K. Yessad: Sep 2010: merge CMAC and CNMI + cleanings.
!      R. El Khatib : 16-Jul-2012 Fullpos move away from STEPO
!      H. Varella: July 2012: Read in ARPEGE format and write in RAW format
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      P. Brousseau (Apr 2014) IAU
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      K. Yessad (July 2014): Move some variables.
!     F. Vana  05-Mar-2015  Support for single precision
!      B. Bochenek (Apr 2015): Phasing: update
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!      A. Geer      27 Jul 2015  More OOPS cleaning: VARBC by argument
!      R. El khatib 22-Feb-2016 NSTEP passed to opdis
!      O. Marsden    Aug 2016  Removed use of SPA3
!      P. Lean       22 Mar 2017 OOPS cleaning: Jo-table by argument
!      T. Montmerle  Jul 2018: Jo-table by argument for CNT3_GLO
!     ------------------------------------------------------------------

USE TYPE_MODEL    , ONLY : MODEL
USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE FIELDS_MOD    , ONLY : FIELDS
USE MTRAJ_MOD     , ONLY : MTRAJ
USE PARKIND1      , ONLY : JPRD, JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN        , ONLY : NULOUT
USE JO_TABLE_MOD  , ONLY : JO_TABLE
USE YOMCT0        , ONLY : NCONF, LELAM
USE YOMCT3        , ONLY : NSTEP
USE YOMTIM        , ONLY : RSTART, RVSTART, RTIMEF
USE GFL_SUBS_MOD  , ONLY : DEACT_CLOUD_GFL, REACT_CLOUD_GFL
USE YOMFP_SERV    , ONLY : FP_SERV_C001
USE VARBC_CLASS   , ONLY : CLASS_VARBC
USE TOVSCV_MOD    , ONLY : TOVSCV
USE SUPERGOM_CLASS, ONLY : CLASS_SUPERGOM
USE DBASE_MOD     , ONLY : DBASE
#if defined(WITH_OASIS) || defined(WITH_NEMO)
USE COUPLING , ONLY : CPL_CONFIG
#endif
USE FULLPOS       , ONLY : TFPOS
USE YOMARG        , ONLY : NUDATE
USE YOMFPC        , ONLY : LFPMOIS

!      -----------------------------------------------------------

IMPLICIT NONE


TYPE(GEOMETRY)      ,INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for call to RERESF
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)         ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
TYPE(JO_TABLE)      ,INTENT(INOUT) :: YDJOT
TYPE(CLASS_VARBC)   ,INTENT(INOUT), OPTIONAL :: YDVARBC
TYPE(TOVSCV)        ,INTENT(IN), OPTIONAL    :: YDTCV
TYPE(CLASS_SUPERGOM),INTENT(INOUT), OPTIONAL :: YDGOM5
CLASS(DBASE)        ,INTENT(INOUT), OPTIONAL :: YDODB
TYPE(TFPOS)      ,INTENT(IN), OPTIONAL :: YDFPOS

INTEGER(KIND=JPIM) :: IMMCLI
LOGICAL :: LLINITMONTH ! .TRUE. to control the climatology month against the month of the initial file

REAL(KIND=JPRD) :: ZCT, ZVT, ZWT

CHARACTER (LEN=9)   :: CLCONF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------

#include "user_clock.intfb.h"

#include "cnt4.intfb.h"
#include "opdis.intfb.h"
#include "su3yom.intfb.h"
#include "cnt3_glo.intfb.h"
#include "cnt3_lam.intfb.h"

#include "fcttim.func.h"

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CNT3',0,ZHOOK_HANDLE)

!      -----------------------------------------------------------
  
!*       1.    Initialize LEVEL 3 COMMONS.
!              ---------------------------

WRITE(UNIT=NULOUT,FMT='('' START CNT3'')')
CALL SU3YOM
IF (YDMODEL%YRML_PHY_EC%YREPHY%LEPHYS) THEN
 !IF (YDMODEL%YRML_PHY_EC%YREPHY%LEPCLD.OR.YDMODEL%YRML_PHY_SLIN%YRPHNC%LEPCLD2) THEN
  IF (YDMODEL%YRML_PHY_EC%YREPHY%LEPCLD) THEN
    CALL REACT_CLOUD_GFL(YDGEOMETRY%YRDIMV,YDGEOMETRY%YRCVER,YDMODEL)
  ELSE
    CALL DEACT_CLOUD_GFL(YDGEOMETRY%YRDIMV,YDGEOMETRY%YRCVER,YDMODEL)
  ENDIF
ENDIF

! * Coupling initialization for OASIS or single executable:
! * NB: This needs to be done before we read the restarts 
#if defined(WITH_OASIS) || defined(WITH_NEMO)
#ifdef WITH_NEMO
YDMODEL%YRML_AOC%YRMCC%CPLNG_ACTIVE=YDMODEL%YRML_AOC%YRMCC%LNEMOCOUP
#endif
CALL CPL_CONFIG(YDGEOMETRY,YDMODEL%YRML_AOC%YRMCC) 
#endif

!      -----------------------------------------------------------

!*       2.    RESTART.
!              --------

IF(NCONF == 1.OR.NCONF == 302) THEN
!  ** OPDIS IS CALLED HERE ONLY TO AVOID MEMORY FRAGMENTAION **
  CLCONF='000000000'
  CALL USER_CLOCK(PELAPSED_TIME=ZWT,PVECTOR_CP=ZVT,PTOTAL_CP=ZCT)
  ZCT=ZCT-RSTART
  ZVT=ZVT-RVSTART
  ZWT=ZWT-RTIMEF
  CALL OPDIS(CLCONF,'CNT3',ZCT,ZVT,ZWT,RSTART,RTIMEF,NSTEP,YDMODEL%YRML_DYN%YRDYNA%LNHDYN)
ENDIF

IF (.NOT.FP_SERV_C001%LFP_SERVER) THEN
  IF (PRESENT(YDFPOS)) THEN
    IF (LFPMOIS) THEN
    ! Date of the starting file (no need to change the climatology file)
      LLINITMONTH=.TRUE.
      IMMCLI=NMM(NUDATE)
    ELSE
      LLINITMONTH=.FALSE.
    ENDIF
  ELSE
    LLINITMONTH=.FALSE.
  ENDIF
  IF (LLINITMONTH) THEN
    IF (LELAM) THEN
      CALL CNT3_LAM(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDVARBC,KINITMONTH=IMMCLI)
    ELSE
      CALL CNT3_GLO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,KINITMONTH=IMMCLI)
    ENDIF
  ELSE
    IF (LELAM) THEN
      CALL CNT3_LAM(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDVARBC)
    ELSE
      CALL CNT3_GLO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDGOM5=YDGOM5,YDODB=YDODB)
    ENDIF
  ENDIF
ENDIF

CALL CNT4(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDGOM5=YDGOM5, &
        & YDODB=YDODB,YDFPOS=YDFPOS)

WRITE(UNIT=NULOUT,FMT='('' END CNT3'')')

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CNT3',1,ZHOOK_HANDLE)
END SUBROUTINE CNT3
