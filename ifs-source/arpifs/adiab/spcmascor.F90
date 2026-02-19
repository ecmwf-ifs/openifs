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

SUBROUTINE SPCMASCOR(YDLAP,YDDYN,PSPSP)

!**** *SPCMASCOR* - MASS CORRECTION IN SPECTRAL SPACE.
!                   "MASCOR" cheap formulation of mass corrector.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCMASCOR(..)

!        Explicit arguments :
!        -------------------- 

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. Yessad (after old part 0.5 of SPCSI, with cleanings)
!        Original : 09-Dec-2004

!     Modifications.
!     --------------
!        simplify : Nils Wedi + Mats Hamrud 2008-02-08
!        K. Yessad (Sep 2008): update comments + cleanings.
!        K. Yessad (Feb 2012): tests in the caller.
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULERR
USE YOMMP0             , ONLY : NPRINTLEV
USE YOMCT3             , ONLY : NSTEP
USE YOMDYN             , ONLY : TDYN
USE YOMLAP             , ONLY : TLAP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TLAP)          ,INTENT(IN)    :: YDLAP
TYPE(TDYN)          ,INTENT(IN)    :: YDDYN
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPSP(:) 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCMASCOR',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    GLOBAL MASS CORRECTION.
!              -----------------------

IF (NPRINTLEV >= 1 ) THEN
  WRITE(NULERR,*) ' '
  WRITE(NULERR,*) '*** SPCMASCOR at step ',NSTEP,' *** START ***'
  WRITE(NULERR,*) '   Initial mass = ',YDDYN%GMASSI
  WRITE(NULERR,*) '   Current mass = ',YDDYN%GMASS0
  WRITE(NULERR,*) '   Increment    = ',YDDYN%GMASSINC
  WRITE(NULERR,*) '   Uncorrected value:'
  WRITE(NULERR,*) '   ln(ps) (m=0,n=0) = ',PSPSP(YDLAP%NASM0(0))
  WRITE(NULERR,*) '   Equiv.  mass = ',EXP(PSPSP(YDLAP%NASM0(0)))
ENDIF

! * apply mass correction.
PSPSP(YDLAP%NASM0(0))=PSPSP(YDLAP%NASM0(0))-YDDYN%GMASSINC

IF ( NPRINTLEV >= 1 ) THEN
  WRITE(NULERR,*) '   Corrected value:'
  WRITE(NULERR,*) '   ln(ps) (m=0,n=0) = ',PSPSP(YDLAP%NASM0(0))
  WRITE(NULERR,*) '   Equiv.  mass = ',EXP(PSPSP(YDLAP%NASM0(0)))
  WRITE(NULERR,*)'*** SPCMASCOR at step ',NSTEP,' ***  END  ***'
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPCMASCOR',1,ZHOOK_HANDLE)
END SUBROUTINE SPCMASCOR
