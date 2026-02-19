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

SUBROUTINE SUMP0(KULOUT,KULNAM)

!**** *SUMP0*   - Sets up distributed memory variables - part 0

!     Purpose.
!     --------
!       * Define basic distributed memory variables:
!         YOMMP0 variables which are in NAMPAR1

!**   Interface.
!     ----------
!        *CALL* *SUMP0*

!        explicit arguments :
!        --------------------

!        implicit arguments :
!        --------------------

!     method.
!     -------
!        see documentation

!     externals.
!     ----------

!     reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     author.
!     -------
!      MPP Group *ECMWF*
!      original : 95-10-01

!     modifications.
!     --------------
!      R. El Khatib  11-May-2007 LSLONDEM=.T. as default
!      Y. Seity      31-Aug-2007 LSYNC_SLCOM LSYNC_TRANS
!      K. Yessad (Jan 2010): remove useless variables.
!      P. Marguinaud 18-May-2010 Add NDISTIO 
!      P. Marguinaud 14-Jun-2010 Add a message about NDISTIO 
!      P. Marguinaud 11-Sep-2012 IO parameters initialization + cleaning
!      P. Marguinaud 10-Oct-2013 Remove various IO options
!      K. Yessad (July 2013): minor modifications in printings (printing and setup in the same routine).
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib  15-Sep-2016 Setup of NOUTTYPE
!      K. Yessad (Dec 2016): Prune obsolete options.
!      R. El Khatib 22-May-2019 LGPTOT_CAP
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM , JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : NCONF, LFDBOP, LALLOPR, LECMWF, LELAM, LARPEGEF, L_OOPS
USE YOMMP0   , ONLY : NPROC, LSPLIT, LSPLITOUT, NWRTOUT, NSTRIN, NSTROUT, &
 & NFLDIN, NSLPAD, NOUTTYPE, LEQ_REGIONS, NCOMBFLEN, LGPTOT_CAP,LVECADIN, &
 & LSYNC_SLCOM, LSYNC_TRANS, M_BARRINC_DIWRGRID, L_GATHERV_WRGP, LSLDEBUG, &
 & NDISTIO, LUSEWRGRIDALL, NPRINTLEV, NPRTRW, MYSETW, NPRTRV, MYSETV,&
 & LSLONDEM,NTRANS_SYNC_LEVEL
USE YOM_GRIB_CODES  , ONLY : NGRBVO, NGRBD, NGRBT, NGRBQ, NGRBO3, NGRBLNSP, NGRBZ, &
 & NGRBCLWC, NGRBCIWC, NGRB118, NGRB119
USE YOMLUN   , ONLY : NULOUT, NULERR
USE SPECTRAL_FIELDS_MOD, ONLY : SETUP_SPEC

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULNAM 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISTAT, IU
LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------

#include "nampar1.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMP0',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!     1. DEFAULT VALUES.
!     ------------------

IF(NCONF/100 == 1) THEN
  LSPLIT=.FALSE.
ELSE
  LSPLIT=.TRUE.
ENDIF
LSLONDEM=.TRUE.
LSPLITOUT=.FALSE.
IF ( LECMWF ) THEN
  LEQ_REGIONS=.TRUE.
ELSE
  LEQ_REGIONS=.FALSE.
ENDIF
IF(LFDBOP) THEN
  NWRTOUT=1
ELSE
  NWRTOUT=NPROC
ENDIF
NSLPAD=0
IF (LARPEGEF) THEN
  NOUTTYPE=1
ELSE
  NOUTTYPE=-1
ENDIF
NFLDIN=0
IF ( LECMWF ) THEN
  NSTRIN=1
  NSTROUT=0
ELSE
  NSTRIN=NPROC
  NSTROUT=1
ENDIF

! set ncombflen (in words) to a value less than the mailbox size (in bytes).
NCOMBFLEN=1800000

LSYNC_SLCOM=.FALSE.
LSYNC_TRANS=.FALSE.
LSLDEBUG=.FALSE.

NDISTIO=0
LUSEWRGRIDALL=.FALSE.

LGPTOT_CAP=LELAM
!     ------------------------------------------------------------------

!     2. Read NAMPAR1 namelist.
!     -------------------------

CALL POSNAM(KULNAM,'NAMPAR1')
READ(KULNAM,NAMPAR1)

!     ------------------------------------------------------------------

!     3. Checkings and resettings.
!     ----------------------------

NSTRIN = MIN (NSTRIN, NPROC)
IF (NSTRIN <=0) NSTRIN = NPROC
NSTROUT = MIN (NSTROUT, NPROC)
IF (NSTROUT <=0) NSTROUT = NPROC

IF (.NOT. (NDISTIO(1) >= 0 .AND. NDISTIO(1) <=3)) THEN
  WRITE(KULOUT,*) " SUMP0: NDISTIO(1) must be either 0, 1, 2 or 3"
  CALL FLUSH(KULOUT)
  CALL ABOR1(' SUMP0: ABOR1 CALLED')
ENDIF

WRITE(KULOUT, *) " NDISTIO(1) = ", NDISTIO(1)
IF (NDISTIO(1) == 1) THEN
  WRITE (KULOUT, '(I6," FILES WILL BE CREATED ON OUTPUT")') NSTROUT
  WRITE (KULOUT, *)
ENDIF

IF( NOUTTYPE > 0) THEN
  IF(NOUTTYPE == 2) THEN
    LFDBOP=.TRUE.
  ELSE
    LFDBOP=.FALSE.
  ENDIF
ENDIF

IF(LFDBOP) THEN
   NOUTTYPE=2
   LSPLITOUT=.TRUE.
ELSE
   NOUTTYPE=1
ENDIF

IF(.NOT.LSPLITOUT) NWRTOUT=NPROC
IF( NSLPAD < 0 ) CALL ABOR1(' SUMP0: NSLPAD must be >=0 ')

!  Setup distributed spectral field type

CALL SETUP_SPEC(NPRTRV, NPRTRW, MYSETV, MYSETW, LELAM, NULOUT, &
 & NGRBVO, NGRBD, NGRBT, NGRBQ, NGRBO3, NGRBLNSP, NGRBZ, &
 & NGRBCLWC, NGRBCIWC, NGRB118, NGRB119, -99)

IF(L_OOPS .AND. LVECADIN) THEN
  ! NVSEPL and NVSEPC not calculated in OOPS
  CALL ABOR1 (' SETUP_VAR: L_OOPS .AND. LVECADIN incompatible')
ENDIF

!     ------------------------------------------------------------------

!     4. Printings.
!     -------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = KULOUT
WRITE(IU,*) ''
WRITE(IU,'('' ----- PRINTINGS IN SUMP0: '')')
WRITE(IU,'('' NDISTIO = '',30I2)') NDISTIO
WRITE(IU,'('' NCOMBFLEN = '',I12)') NCOMBFLEN
WRITE(IU,'('' LSYNC_SLCOM = '',L2)') LSYNC_SLCOM
WRITE(IU,'('' LSYNC_TRANS = '',L2)') LSYNC_TRANS
WRITE(IU,'('' LSLONDEM = '',L2)')LSLONDEM
WRITE(IU,'('' NSTRIN = '',I4,'' NSTROUT = '',I4)') NSTRIN,NSTROUT  
WRITE(IU,'('' NOUTTYPE = '',I2)') NOUTTYPE
WRITE(IU,'('' NFLDIN = '',I4)') NFLDIN
WRITE(IU,'('' LSPLITOUT = '',L2)') LSPLITOUT
WRITE(IU,*) ' NWRTOUT = ',NWRTOUT
WRITE(IU,*) ' M_BARRINC_DIWRGIRD = ',M_BARRINC_DIWRGRID
WRITE(IU,'('' L_GATHERV_WRGP = '',L2)') L_GATHERV_WRGP
WRITE(IU,'('' LUSEWRGRIDALL = '',L2)') LUSEWRGRIDALL
WRITE(IU,'('' LEQ_REGIONS = '',L2)') LEQ_REGIONS
WRITE(IU,'('' LSPLIT = '',L2)') LSPLIT
WRITE(IU,'('' LSLDEBUG = '',L2)') LSLDEBUG
WRITE(IU,'(''  LVECADIN = '',L2)') LVECADIN
IF (LELAM) THEN
  WRITE(IU,'('' LGPTOT_CAP = '',L2)') LGPTOT_CAP
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMP0',1,ZHOOK_HANDLE)
END SUBROUTINE SUMP0
