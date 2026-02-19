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

SUBROUTINE SULUN

!**** *SULUN * - Routine to initialize the common YOMLUN

!     Purpose.
!     --------
!           Initialize and print the common YOMLUN

!**   Interface.
!     ----------
!        *CALL* *SULUN

!        Explicit arguments :
!        --------------------

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      K. Yessad (Aug 2009): prune conf 912.
!      P. Marguinaud 01-Jan-2011 IO server LUN reservation.
!      P. Marguinaud 10-Oct-2013 Cleaning
!      K. Yessad (Oct 2013): call to SUMPOUT moved in SUMPINI_PRT.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (Dec 2016): Prune obsolete options.
!      R. El Khatib : 02-Aug-2018 argument sumpout
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY :  NULOUT   ,NULNAM   ,NPOSSH   ,&
 & NTIDE    ,NTRJSH   ,NINISH   ,NINIGG   ,NFGISH   ,&
 & NPPPSH   ,NPODDH   ,NULCL1   ,NULCL2   ,&
 & NULASE   ,NULASS   ,NULDILA  ,NULCONT  ,&
 & NULRCF   ,NULHWF   ,NULUSR1  ,NULUSR2  ,&
 & NULUSR3  ,NULUSR4  ,NULUSR5  ,NULCO    ,NEFLS    ,&
 & NEFLSS   ,&
 & NULFP01  ,NULFP02  ,NULFP03  ,NULFP04  ,NULFP05  ,&
 & NULFP06  ,NULFP07  ,NULFP08  ,NULFP09  ,NULFP10  ,&
 & NULFP11  ,NULFP12  ,NULFP13  ,NULFP14  ,NULFP15  ,&
 & NULFPOS  ,NSCRTCH  ,NULERR   ,&
 & NULRAD   ,NULRTL   ,&
 & NEGASH   ,NULTRAJHR,NULTRAJBG,&
 & NUIO_SERV_LOG
USE YOMMP0   , ONLY : LSCMEC

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "sumpout.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SULUN',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!        1.    Initialize YOMLUN.
!              ------------------

NULOUT = 20 ! Note: NULOUT was set to 6 in YOMLUN_IFSAUX. Reset here.

!     If DM version diagnostic output optionally only on subset of PEs
! KULOUT.XML => LUN=40, OPENS IN GSTAT !

IF (.NOT. LSCMEC) THEN
  CALL SUMPOUT('NODE')
ENDIF

!     ------------------------------------------------------------------

!        2.    Print YOMLUN.
!              -------------

WRITE(UNIT=NULOUT,FMT='(&
 & '' NULOUT='',I3,'' NULNAM='',I3,/,&
 & '' NPOSSH='',I3,'' NPODDH='',I3,'' NULCO='',I3,&
 & '' NTIDE ='',I3  )')&
 & NULOUT,NULNAM &
 & ,NPOSSH,NPODDH,NULCO &
 & ,NTIDE  
WRITE(UNIT=NULOUT,FMT='(&
 & '' NULDILA='',I3,'' NULCONT='',I3,&
 & '' NINISH='',I3,'' NINIGG='',I3,&
 & '' NFGISH='',I3,/,&
 & '' NTRJSH='',I3,&
 & '' NPPPSH='',I3,&
 & '' NULCL1='',I3,'' NULCL2='',I3,/,&
 & '' NULASE='',I3,'' NULASS='',I3,&
 & '' NULRCF='',I3,'' NULHWF='',I3)')&
 & NULDILA,NULCONT &
 & ,NINISH,NINIGG,NFGISH &
 & ,NTRJSH,NPPPSH,NULCL1,NULCL2 &
 & ,NULASE,NULASS,NULRCF,NULHWF  
WRITE(UNIT=NULOUT,FMT='(&
 & '' NULUSR1 = '',I3,'' NULUSR2 = '',I3,'' NULUSR3 = '',I3,&
 & '' NULUSR4 = '',I3,'' NULUSR5 = '',I3)')&
 & NULUSR1,NULUSR2,NULUSR3,NULUSR4,NULUSR5  
WRITE(UNIT=NULOUT,FMT='(&
 & '' NEFLS = '',I3,'' NEFLSS = '',I3)')&
 & NEFLS,NEFLSS
WRITE(UNIT=NULOUT,FMT='(&
 & '' NULFP01 = '',I3,'' NULFP02 = '',I3,'' NULFP03 = '',I3,&
 & '' NULFP04 = '',I3,'' NULFP05 = '',I3,'' NULFP06 = '',I3)')&
 & NULFP01, NULFP02, NULFP03, NULFP04, NULFP05, NULFP06  
WRITE(UNIT=NULOUT,FMT='(&
 & '' NULFP07 = '',I3,'' NULFP08 = '',I3,'' NULFP09 = '',I3,&
 & '' NULFP10 = '',I3,'' NULFPOS = '',I3,'' NSCRTCH = '',I3)')&
 & NULFP07, NULFP08, NULFP09, NULFP10, NULFPOS, NSCRTCH  
WRITE(UNIT=NULOUT,FMT='(&
 & '' NULERR = '',I3,&
 & '' NULRAD = '',I3,'' NULRTL = '',I3)')&
 & NULERR, NULRAD, NULRTL  
WRITE(UNIT=NULOUT,FMT='(&
 & '' NULFP11 = '',I3,'' NULFP12 = '',I3,'' NULFP13 = '',I3,&
 & '' NULFP14 = '',I3,'' NULFP15 = '',I3)')&
 & NULFP11, NULFP12, NULFP13, NULFP14, NULFP15  
WRITE(UNIT=NULOUT,FMT='('' NEGASH = '',I3)') NEGASH
WRITE(UNIT=NULOUT,FMT='('' NULTRAJHR = '',I3)') NULTRAJHR
WRITE(UNIT=NULOUT,FMT='('' NULTRAJBG = '',I3)') NULTRAJBG
WRITE(UNIT=NULOUT,FMT='('' NUIO_SERV_LOG = '',I3)') NUIO_SERV_LOG

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SULUN',1,ZHOOK_HANDLE)
END SUBROUTINE SULUN
