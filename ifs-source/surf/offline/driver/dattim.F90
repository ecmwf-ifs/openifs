SUBROUTINE DATTIM(PJUL,KYMD,KHM)
USE PARKIND1  ,ONLY : JPIM     ,JPRB , JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
! (C) Copyright 1996- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *DATTIM*  - Gives time of the model in the format yyyymmdd and hhmm

!     Purpose.
!     --------
!     Gives time of the model

!**   Interface.
!     ----------
!        *CALL* *DATTIM

!        Explicit arguments :
!        --------------------

!     KYMD   CURRENT DATE IN THE MODEL (YYYYMMDD)      OUTPUT
!     KHM    CURRENT TIME IN THE MODEL (HHMM)          OUTPUT
!     PJUL   Current fractional Julian day, i.e. a counter in days
!             that corresponds to 0 at 1 January of the first year
!             of simulation

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------

!     Author.
!     -------
!        Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 96-03-03

!     ------------------------------------------------------------------


USE YOMDYN1S , ONLY : NSTEP    ,TSTEP
USE YOMCST   , ONLY : RDAY
USE YOMRIP   , ONLY : RTIMST   ,NINDAT   ,NSSSSS

IMPLICIT NONE

!* arguments
!
REAL(KIND=JPRD),INTENT(OUT)    :: PJUL 
INTEGER(KIND=JPIM),INTENT(OUT)    :: KYMD 
INTEGER(KIND=JPIM),INTENT(OUT)    :: KHM 

!* Local variables
!
REAL(KIND=JPRD) :: ZTIME,ZTIMCUR,ZTIM1JAN
INTEGER(KIND=JPIM) :: ISS,IHH,IMM,IYY,IDAY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttim.h"
#include "incdat.intfb.h"

IF (LHOOK) CALL DR_HOOK('DATTIM',0,ZHOOK_HANDLE)


!     ------------------------------------------------------------------

!*       1.   TIME OF THE MODEL.
!             ------------------

ZTIME=REAL(NSTEP,KIND=JPRD)*TSTEP
ISS=MOD(INT(ZTIME)+NSSSSS,86400)
IHH=ISS/3600
IMM=MOD(ISS,3600)/60
KHM=100*IHH+IMM
! ztimcur is the julian time at the beginning of the time step
ztimcur=RTIMST+ztime

!*       2.   DATE OF THE MODEL.
!             ------------------

IDAY=(INT(ZTIME)+NSSSSS)/86400
CALL INCDAT(NINDAT,IDAY,KYMD)

!CBH!*       3.   Fractional Julian day, with 0 being 1 January 00 UT, in the
!CBH!             first year of the simulation
!*       3.   Fractional Julian day, with 1 being 1 January 00 UT, in the
!             first year of the simulation

iyy=nccaa(nindat)
ztim1jan=rtime(iyy,1,1,0)
!CBH      pjul=(ztimcur-ztim1jan+1.)/rday
pjul=(ztimcur-ztim1jan)/rday+1.

IF (LHOOK) CALL DR_HOOK('DATTIM',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE DATTIM
