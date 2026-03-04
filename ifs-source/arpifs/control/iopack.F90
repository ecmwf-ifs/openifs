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

SUBROUTINE IOPACK(YDGEOMETRY,YDFIELDS,YDGFL5,YDMODEL,CDCONF,PTRAJEC,YDECV)

!**** *IOPACK*  - input/output handling during integration

!     Purpose.
!     --------
!     Write or read the trajectory, write history file ,write out
!     post-processed data or create restart files.

!**   Interface.
!     ----------
!        *CALL* *IOPACK(CDCONF)

!        Explicit arguments :     CDCONF - configuration of call
!        --------------------

!        Implicit arguments :      None.
!        --------------------

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
!      David Dent                         *ECMWF*
!      Original : 88-08-02

!     Modifications.
!     --------------
!      Y.Seity      11-Jan-2008 add surfex output files (LLS)
!      A.Alias      07-Aug-2009 modify LLA/LLR/LLS for surfex output files ('L'/'Q')
!      K. Yessad: Sep 2010 : simplify organigramme.
!      R. El Khatib : 01-Mar-2012 LFPOS => LECFPOS
!      R. El Khatib 17-Jul-2012 Fullpos move away
!      P.Marguinaud Oct 2013 Move SURFEX IO after model IO
!      F. Vana  28-Nov-2013 : Redesigned trajectory handling
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden  Aug 2016 : Remove use of SPA3
!      R. El Khatib : 23-Aug-2016 lgrbop <=> larpegef
!      K. Yessad (Dec 2016): Prune obsolete options.
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE YOMGFL             , ONLY : TGFL
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : LARPEGEF, NCONF, LFDBOP, LELAM, LECFPOS,L_OOPS
USE YOMCT3             , ONLY : NSTEP
USE YOMVRTL            , ONLY : LTLINT
USE YOMMP0             , ONLY : NOUTTYPE
USE YOMSENS            , ONLY : NJROPT
USE YOMVAR             , ONLY : LTRREF
USE YOMTRAJ            , ONLY : TRAJ_TYPE
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE TYPE_ECV           , ONLY : ECV_CONTAINER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for call to RD801
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL5
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF
TYPE(TRAJ_TYPE), OPTIONAL, INTENT(IN)    :: PTRAJEC
TYPE(ECV_CONTAINER),  OPTIONAL, INTENT(IN)    :: YDECV
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IOUTTYPE

LOGICAL :: LL801, LLA, LLB, LLFDBOP,&
 & LLI, LLJ, LLL, LLR, LLS, LLT, LLV
CHARACTER(LEN=9) :: CLIOPACK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "aro_surf_diagh.h"

#include "abor1.intfb.h"
#include "add3to5.intfb.h"




#include "ecoptra.intfb.h"
#include "suspe0.intfb.h"
#include "wrmlpp.intfb.h"
#include "wrmlppa.intfb.h"
#include "wrmlpplg.intfb.h"
#include "wrresf.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('IOPACK',0,ZHOOK_HANDLE)

CLIOPACK='IOPACK<'//CDCONF//'>'
ASSOCIATE(YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)
ASSOCIATE(NFLEVL=>YDGEOMETRY%YRDIMV%NFLEVL, &
 & NSTOP=>YDRIP%NSTOP)
!     ------------------------------------------------------------------

!*       1.    WRITE TRAJECTORY.
!              -----------------

CALL GSTATS(11,0)


! 39r1 still unused letters: C,G,H,N,O,P,X
! 40r1 currently unused letters: C,E,G,H,K,M,N,O,P,U,X,Y,Z
! cy44: currently unused letters: C,E,G,H,K,M,N,O,P,U,W,X,Y,Z
IF (.NOT. ANY(SPREAD(CDCONF,1,13)==(/'A','B','D','F','I','J','L','Q','R','S','T','V','f'/))) THEN
  WRITE(0,*) "IOPACK: CDCONF=",CDCONF
  CALL ABOR1(' IOPACK: wrong value for CDCONF ')
ENDIF

LLA  = CDCONF == 'A'.OR.CDCONF == 'F'.OR.CDCONF == 'D'.OR.CDCONF == 'Q'.OR.CDCONF == 'f'
LLB  = CDCONF == 'B'
LLI  = CDCONF == 'I'
LLJ  = CDCONF == 'J'
LLR  = CDCONF == 'F'.OR.CDCONF == 'R'.OR.CDCONF == 'L'.OR.CDCONF == 'f'
LLS  = CDCONF == 'S'.OR.CDCONF == 'Q'.OR.CDCONF == 'L'.OR.CDCONF == 'f'
LLT  = CDCONF == 'T'
LLL  = CDCONF == 'L'
LLV  = CDCONF == 'V'
LL801=NSTEP == NSTOP.AND.NCONF == 801.AND..NOT.LTLINT

IF(.NOT.L_OOPS .AND. NCONF/100 == 1) THEN
  LLFDBOP=LFDBOP
  LFDBOP=.FALSE.
  IOUTTYPE=NOUTTYPE
  NOUTTYPE=1
ENDIF

IF (LLT) THEN
  CALL COPTRA(CDCONF)
  IF (LELAM) CALL ECOPTRA(YDGEOMETRY%YRDIMV,YDRIP,CDCONF,YDFIELDS%YRSPEC)
ENDIF

!     ------------------------------------------------------------------

!*       2.   READ TRAJECTORY.
!             ----------------

IF (LLT) THEN
! COPTRA doing nothing
  CALL COPTRA(CDCONF)
  IF (LELAM) CALL ECOPTRA(YDGEOMETRY%YRDIMV,YDRIP,CDCONF,YDFIELDS%YRSPEC)
ENDIF

IF(LLB) THEN
  IF (LL801.AND. NJROPT == 1) THEN
    IF (.NOT.LTRREF) THEN
!     COPTRA doing nothing
      CALL COPTRA(CDCONF)
      IF (LELAM) CALL ECOPTRA(YDGEOMETRY%YRDIMV,YDRIP,CDCONF,YDFIELDS%YRSPEC)
    ELSE
      CALL RD801(YDGEOMETRY,YDFIELDS%YRGFL,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYN,YDDYNA,YDMODEL%YRML_LBC, &
 &               YDFIELDS%YRSPEC)
    ENDIF
  ELSE
!   COPTRA doing nothing
    CALL COPTRA(CDCONF)
    IF (LELAM) CALL ECOPTRA(YDGEOMETRY%YRDIMV,YDRIP,CDCONF,YDFIELDS%YRSPEC)
  ENDIF
ENDIF
IF(LLV) THEN
  IF((NCONF == 801).AND.(NJROPT == 1))THEN
    IF(.NOT. LTRREF) THEN
      CALL COPTRA(CDCONF)
      IF (LELAM) THEN
        CALL ECOPTRA(YDGEOMETRY%YRDIMV,YDRIP,CDCONF,YDFIELDS%YRSPEC)
        CALL ADD3TO5(YDGEOMETRY%YRMP,YDMODEL%YRML_GCONF,YDFIELDS%YRSPEC)
        CALL SUSPE0(YDGEOMETRY,YDMODEL%YRML_GCONF%YRDIMF,YDFIELDS%YRSPEC)
        YDFIELDS%YRSPEC%SP2D(1:NFLEVL,1:2)=0.0_JPRB
      ENDIF
    ELSE
      CALL RD801(YDGEOMETRY,YDFIELDS%YRGFL,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN%YRDYN,YDDYNA, &
 &               YDMODEL%YRML_LBC,YDFIELDS%YRSPEC)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    WRITE RESTART FILES
!              -------------------

IF (LLR) THEN
  CALL WRRESF(YDGEOMETRY,YDFIELDS,YDMODEL)
ENDIF

!     ------------------------------------------------------------------

!*       4.    WRITE MODEL SURFACE DATA.
!              -----------------------

IF (LLS) THEN
  CALL ARO_SURF_DIAGH(YDGEOMETRY,YDMODEL,YDFIELDS%YRGFL,YDFIELDS%YRSURF,YDFIELDS%YRSPEC,YDFIELDS%YRCFU,YDFIELDS%YRXFU,YDRIP)
ENDIF

!     ------------------------------------------------------------------

!*       5.    WRITE MODEL LEVEL DATA.
!              -----------------------

IF (LLA) THEN
  IF (LELAM) THEN
    CALL WRMLPPA(YDGEOMETRY,YDFIELDS%YRGFL,YDFIELDS%YRSURF,YDFIELDS%YRSPEC,YDFIELDS%YRCFU,YDFIELDS%YRXFU,YDMODEL,CDCONF,YDFIELDS%YMCUF)
  ELSE
    CALL WRMLPP(YDGEOMETRY,YDFIELDS%YRGFL,YDGFL5,YDFIELDS%YRSURF,YDFIELDS%YRSPEC,YDFIELDS%YRCFU,YDFIELDS%YRXFU, &
 &              YDMODEL,CDCONF,PTRAJEC=PTRAJEC,YDMCUF=YDFIELDS%YMCUF,YDECV=YDECV)
  ENDIF
ENDIF

IF (LLL) THEN
  IF (.NOT.LARPEGEF) THEN
    IF (LECFPOS) THEN
!     CALL WRMLFPL(CDCONF)
      CALL ABOR1('IOPACK : WRMLFPL NO LONGER SUPPORTED')
    ELSE
      CALL WRMLPPLG(YDGEOMETRY,YDFIELDS%YRSURF,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_PHY_SLIN%YREPHLI,YDRIP, &
                  & YDMODEL%YRML_DYN%YRDYN,YDMODEL%YRML_PHY_MF%YRPHY)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       8.    WRITE INCREMENTS
!              -------------------

!IF (LLI) THEN
!  CALL ABOR1('IOPACK : CONFIGURATION ''I'' NOT EXPECTED')
!ENDIF

!     ------------------------------------------------------------------

!*       9.    READ INCREMENTS.
!              ----------------

IF (LLJ) THEN
  CALL ABOR1('IOPACK : CONFIGURATION ''J'' NOT EXPECTED')
ENDIF

!     ------------------------------------------------------------------

IF(NCONF/100 == 1) THEN
  LFDBOP=LLFDBOP
  NOUTTYPE=IOUTTYPE
ENDIF

CALL GSTATS(11,1)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK(CLIOPACK,1,ZHOOK_HANDLE)
END SUBROUTINE IOPACK
