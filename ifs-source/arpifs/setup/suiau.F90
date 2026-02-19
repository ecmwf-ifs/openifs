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

SUBROUTINE SUIAU(YDRIP)

!**** *SUIAU*  - Initialize IAU handling, plus basic tests

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUIAU

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See purpose above

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------

!     Author.
!     -------
!        Pierre BROUSSEAU *CNRM/GMAP*
!        Original :    fev-2014

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   ,NULNAM
USE YOMIAU   , ONLY : LIAU     ,ALPHAIAU,NSTARTIAU, NSTOPIAU,TSTARTIAU,TSTOPIAU
USE YOMCT0   , ONLY : NCONF    ,LELAM
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(TRIP),INTENT(INOUT):: YDRIP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namiau.nam.h"
#include "abor1.intfb.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUIAU',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------

WRITE(NULOUT,*) ' *** SUIAU : SETUP IAU Term LAM forecast'

!     1. DEFAULT VALUES
!     -----------------

LIAU=.FALSE.
ALPHAIAU=1.0_JPRB
TSTARTIAU = 0.0_JPRB
TSTOPIAU = 1.0_JPRB
NSTARTIAU = 0_JPIM
NSTOPIAU = 1_JPIM

!     2. READ NAMELIST AND CHECK CONSISTENCY
!     --------------------------------------

CALL POSNAM(NULNAM,'NAMIAU')
READ(NULNAM,NAMIAU)

IF (LIAU) THEN
  IF (.NOT.LELAM) THEN
     CALL ABOR1('IAU not yet coded for global model ')
  ENDIF
  IF (NCONF /= 1) THEN
     CALL ABOR1('LIAU=.TRUE. only with NCONF==1 ')
  ENDIF
  IF (TSTARTIAU > TSTOPIAU) THEN
        CALL ABOR1('TSTARTIAU>TSTOPIAU, CHECK IAU NAMELIST')
  ENDIF

  NSTARTIAU=NINT(TSTARTIAU/TSTEP)
  NSTOPIAU=NINT(TSTOPIAU/TSTEP)

  WRITE(NULOUT,*) 'LIAU switched on'
  WRITE(NULOUT,*) 'part of the increment read in file added during the run :',ALPHAIAU
  WRITE(NULOUT,*) 'between time ',TSTARTIAU,' (sec.) corresponding to timestep ',NSTARTIAU
  WRITE(NULOUT,*) 'and time ',TSTOPIAU,' (sec.) corresponding to timestep ',NSTOPIAU
  WRITE(NULOUT,*) 'part of the increment at each time step :1/',(NSTOPIAU-NSTARTIAU+1)/ALPHAIAU
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUIAU',1,ZHOOK_HANDLE)
END SUBROUTINE SUIAU
