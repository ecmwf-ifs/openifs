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
SUBROUTINE SURIP(YDDIM,YDDYNA,YDRIP,PTSTEP)

!**** *SURIP * - Routine to initialize the module YOMRIP.

!     Purpose.
!     --------
!       Routine to initialize the module YOMRIP (and TSTEP_TRAJ in YOMTRAJ).
!       Timestep, number of timesteps, and some timestep dependent variables.
!       In particular, all variables updated in UPDTIM must be setup there.
!       For OOPS, all calculations are model-dependent ones.

!**   Interface.
!     ----------
!        *CALL* *SURIP(...)

!     Explicit arguments :
!     --------------------
!        none

!     Implicit arguments :
!     --------------------

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
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified : 2011-Mar A.Alias LASTRF added to prevent drift in insolation
!        Y. Bouteloup (Feb 2011) : Add RCODECF  ,RSIDECF  ,RCOVSRF  ,RSIVSRF
!      K. Yessad (July 2014): Model-independent variables moved in YOMRIP0/SURIP0, encapsulation.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      O. Marsden      Nov 2017 Removal of NUSTOP and UTSTEP from YOMARG, replaced by CSTOP and TSTEP in YOMRIP
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM, NULOUT, NULERR
USE YOMCT0   , ONLY : NCONF, L4DVAR
USE YOMRIP0  , ONLY : RTIMST
USE YOMRIP   , ONLY : TRIP
USE YOMTRAJ  , ONLY : TSTEP_TRAJ
USE YOMDYNA  , ONLY : TDYNA
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN)            :: YDDIM
TYPE(TDYNA), INTENT(IN)            :: YDDYNA
TYPE(TRIP) , INTENT(INOUT), TARGET :: YDRIP
REAL(KIND=JPRB), INTENT(IN),OPTIONAL :: PTSTEP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB)   , POINTER :: TSTEP
INTEGER(KIND=JPIM), POINTER :: NSTOP
INTEGER(KIND=JPIM), POINTER :: NFOST
CHARACTER(LEN=8),   POINTER :: CSTOP
INTEGER(KIND=JPIM) :: ISECSPERDAY = 3600*24, ISECSPERHOUR = 3600

#include "namrip.nam.h"

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "sudefo_tstep.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURIP',0,ZHOOK_HANDLE)
ASSOCIATE(NSTADD=>YDRIP%NSTADD, NSTART=>YDRIP%NSTART, NSTASS=>YDRIP%NSTASS, &
 & RCODEC=>YDRIP%RCODEC, RCODECF=>YDRIP%RCODECF, RCODECLU=>YDRIP%RCODECLU, &
 & RCOVSR=>YDRIP%RCOVSR, RCOVSRF=>YDRIP%RCOVSRF, RCOVSRLU=>YDRIP%RCOVSRLU, &
 & RDEASO=>YDRIP%RDEASO, RDECLI=>YDRIP%RDECLI, RDECLU=>YDRIP%RDECLU, &
 & RDTS22=>YDRIP%RDTS22, RDTS62=>YDRIP%RDTS62, RDTSA=>YDRIP%RDTSA, &
 & RDTSA2=>YDRIP%RDTSA2, REQTIM=>YDRIP%REQTIM, RHGMT=>YDRIP%RHGMT, &
 & RIP0=>YDRIP%RIP0, RIP0LU=>YDRIP%RIP0LU, RSIDEC=>YDRIP%RSIDEC, &
 & RSIDECF=>YDRIP%RSIDECF, RSIDECLU=>YDRIP%RSIDECLU, RSIVSR=>YDRIP%RSIVSR, &
 & RSIVSRF=>YDRIP%RSIVSRF, RSIVSRLU=>YDRIP%RSIVSRLU, RSOVR=>YDRIP%RSOVR, & 
 & RSTATI=>YDRIP%RSTATI, RTDT=>YDRIP%RTDT, RTIMTR=>YDRIP%RTIMTR, &
 & RTMOLT=>YDRIP%RTMOLT, RWSOVR=>YDRIP%RWSOVR, TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

! Associate pointers for variables in namelist
TSTEP => YDRIP%TSTEP
NSTOP => YDRIP%NSTOP
NFOST => YDRIP%NFOST
CSTOP => YDRIP%CSTOP
!     ------------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------
! Timesteps
TSTEP=600.0_JPRB
CALL SUDEFO_TSTEP(YDDIM,YDDYNA,YDRIP)
TSTEP_TRAJ=-1._JPRB  ! TSTEP can still be modified

! First timestep
NSTART=0

!! default value
CSTOP = '-9' 

NFOST=0

!     ------------------------------------------------------------------

!*       2.    READ NAMELIST.
!              --------------

CALL POSNAM(NULNAM,'NAMRIP')
READ(NULNAM,NAMRIP)

IF(PRESENT(PTSTEP)) THEN
  IF(PTSTEP /= TSTEP) THEN
    WRITE(NULERR,*) 'SURIP:PTSTEP=',PTSTEP,' TSTEP=',TSTEP
    CALL ABOR1('SURIP:TSTEP FROM NAMELIST DIFFERENT FROM OOPS/JSON')
  ENDIF
ENDIF
!*    2.2
!*    Compute NSTOP from CSTOP and TSTEP 
!    -------------------------------------------
IF (CSTOP /= '-9') THEN
  IF(CSTOP(1:1) == 't') THEN
    READ(UNIT=CSTOP(2:),FMT='(I5)') NSTOP
  ELSE
    IF(TSTEP > 0.0_JPRB) THEN
      READ(UNIT=CSTOP(2:),FMT='(I5)') NSTOP
      IF(CSTOP(1:1) == 'h') THEN
        NSTOP=NINT(NSTOP*ISECSPERHOUR/TSTEP)
      ELSEIF(CSTOP(1:1) == 'd') THEN
        NSTOP=NINT(NSTOP*ISECSPERDAY/TSTEP)
      ELSE
        WRITE(NULOUT,*) ' WRONG FORMAT FOR THE FORECAST'
        WRITE(NULOUT,*) ' SPAN.'
        WRITE(NULOUT,*) ' IT MUST START WITH EITHER'
        WRITE(NULOUT,*) '     t  FOR TIME STEPS or'
        WRITE(NULOUT,*) '     d  FOR DAYS or'
        WRITE(NULOUT,*) '     h  FOR HOURS'
        CALL ABOR1('SURIP: ABOR1 CALLED')
      ENDIF
    ELSE
      WRITE(NULOUT,*) '  ERROR'
      WRITE(NULOUT,*) '  IF YOU SUPPLY THE FORECAST'
      WRITE(NULOUT,*) '  SPAN OTHER THAN IN TIME STEPS'
      WRITE(NULOUT,*) '  YOU HAVE TO SUPPLY THE VALUE'
      WRITE(NULOUT,*) '  OF THE TIME STEP WITH -t ...'
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3     CHECKINGS AND ADDITIONAL SETUP
!              ------------------------------

!*       3.1  Timesteps.

IF (TSTEP_TRAJ<0) TSTEP_TRAJ=TSTEP
IF (TSTEP /= 0._JPRB) THEN
  IF (MOD(TSTEP_TRAJ,TSTEP)/=0) CALL ABOR1('SURIP: TSTEP_TRAJ has to be a multiple of TSTEP')
ENDIF

TDT=TSTEP

!*       3.2  Number of timesteps.

IF (L4DVAR .AND. .NOT.(NCONF==131 .AND. NSTOP>1)) THEN
  ! L4DVAR can be T only if (NCONF=131 and NSTOP>1).
  CALL ABOR1('SURIP: Inconsistency between L4DVAR and NSTOP!')
ENDIF

!*       3.3  Setup variables which are updated by UPDTIM.

NSTADD=0
NSTASS=0
RSTATI=0._JPRB
RTIMTR=RTIMST
RHGMT=0._JPRB
REQTIM=0._JPRB
RSOVR=0._JPRB
RDEASO=0._JPRB
RDECLI=0._JPRB
RWSOVR=0._JPRB
RIP0=0._JPRB
RCODEC=0._JPRB
RSIDEC=0._JPRB
RCOVSR=0._JPRB
RSIVSR=0._JPRB
RCODECF=0._JPRB
RSIDECF=0._JPRB
RCOVSRF=0._JPRB
RSIVSRF=0._JPRB

RDECLU=0._JPRB
RTMOLT=0._JPRB
RIP0LU=0._JPRB
RCODECLU=0._JPRB
RSIDECLU=0._JPRB
RCOVSRLU=0._JPRB
RSIVSRLU=0._JPRB

RDTSA=0._JPRB
RDTSA2=0._JPRB
RDTS62=0._JPRB
RDTS22=0._JPRB

RTDT=0._JPRB

!      ----------------------------------------------------------------

!*       4.    PRINT NAMELIST VARIABLES.
!              -------------------------

WRITE(UNIT=NULOUT,FMT='('' Printings for module YOMRIP '')')
WRITE(UNIT=NULOUT,FMT='('' TSTEP  = '',E14.8,'' TDT    = '',E14.8)') TSTEP,TDT
WRITE(UNIT=NULOUT,FMT='('' TSTEP_TRAJ = '',E14.8)') TSTEP_TRAJ
WRITE(UNIT=NULOUT,FMT='('' NSTART = '',I6)') NSTART
WRITE(UNIT=NULOUT,FMT='('' NSTOP  = '',I6)') NSTOP
WRITE(UNIT=NULOUT,FMT='('' NFOST  = '',I6)') NFOST

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURIP',1,ZHOOK_HANDLE)
END SUBROUTINE SURIP
