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

SUBROUTINE SUDIM_TRAJ(YDDIM)

!------------------------------------------------------------------------------
!**** *SUDIM_TRAJ*   - Initialize geometry (dimensions) for trajectory and background

!     Purpose.
!     --------
!           Initialize some dimensions in YOMTRAJ, for trajectory and background.

!**   Interface.
!     ----------
!        *CALL* *SUDIM_TRAJ

!        Explicit arguments :
!        --------------------
!         none

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

!     Author.
!     -------
!      K. Yessad, after some code previously in SUDYN.
!      Original : July 2014

! Modifications
! -------------
!      E. Holm 22 Jul 2016 : NSMAX_BACKGR04-10 hacked, better later.
! End Modifications
!-------------------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMLUN, ONLY: NULOUT, NULNAM
!! USE TRAJECTORY_MOD, ONLY: NSMAX_TRAJ, NSMAX_BACKGR00, NSMAX_BACKGR01, NSMAX_BACKGR02, NSMAX_BACKGR03
USE YOMTRAJ, ONLY: NSMAX_TRAJ, NSMAX_BACKGR00, NSMAX_BACKGR01, NSMAX_BACKGR02, NSMAX_BACKGR03, &
 & NSMAX_BACKGR04, NSMAX_BACKGR05, NSMAX_BACKGR06, NSMAX_BACKGR07, NSMAX_BACKGR08, NSMAX_BACKGR09, NSMAX_BACKGR10

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN) :: YDDIM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namdim_traj.nam.h"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDIM_TRAJ',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX)
!     ------------------------------------------------------------------


!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

NSMAX_TRAJ=NSMAX
NSMAX_BACKGR00=NSMAX
NSMAX_BACKGR01=NSMAX
NSMAX_BACKGR02=NSMAX
NSMAX_BACKGR03=NSMAX
NSMAX_BACKGR04=NSMAX
NSMAX_BACKGR05=NSMAX
NSMAX_BACKGR06=NSMAX
NSMAX_BACKGR07=NSMAX
NSMAX_BACKGR08=NSMAX
NSMAX_BACKGR09=NSMAX
NSMAX_BACKGR10=NSMAX

!     ------------------------------------------------------------------

!*       2.    Modify default values, read namelist.
!              -------------------------------------

CALL POSNAM(NULNAM,'NAMDIM_TRAJ')
READ(NULNAM,NAMDIM_TRAJ)

!     ------------------------------------------------------------------

!*       3.    Reset variables and test.
!              -------------------------

IF (NSMAX_TRAJ>NSMAX) CALL ABOR1('SUDIM_TRAJ: nsmax_traj>nsmax not possible')

!     ------------------------------------------------------------------

!*       4.    Printings.
!              ----------

WRITE(NULOUT,*)''
WRITE(NULOUT,*)' ===== Printings in SUDIM_TRAJ: '
WRITE(NULOUT,FMT='(''  NSMAX_TRAJ = '',I7)') NSMAX_TRAJ
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR00 = '',I7)') NSMAX_BACKGR00
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR01 = '',I7)') NSMAX_BACKGR01
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR02 = '',I7)') NSMAX_BACKGR02
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR03 = '',I7)') NSMAX_BACKGR03
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR04 = '',I7)') NSMAX_BACKGR04
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR05 = '',I7)') NSMAX_BACKGR05
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR06 = '',I7)') NSMAX_BACKGR06
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR07 = '',I7)') NSMAX_BACKGR07
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR08 = '',I7)') NSMAX_BACKGR08
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR09 = '',I7)') NSMAX_BACKGR09
WRITE(NULOUT,FMT='(''  NSMAX_BACKGR10 = '',I7)') NSMAX_BACKGR10

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDIM_TRAJ',1,ZHOOK_HANDLE)
END SUBROUTINE SUDIM_TRAJ
