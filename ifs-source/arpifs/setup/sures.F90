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

SUBROUTINE SURES(YDRIP,KULOUT)

!**** *SURES*   - Routine to initialize restart control

!     Purpose.
!     --------
!           Initialize restart control common
!**   Interface.
!     ----------
!        *CALL* *SURES(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMRES

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
!      David Dent     *ECMWF*
!      Original : 92-06-01

!     Modifications.
!     --------------
!      Modified by D.Dent       : 02-06-10 LFASTRES option
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified by R. El Khatib : 04-08-12 : print LFASTRES
!      Modified by R. El Khatib : 04-10-12 : path for cp in movies restart
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 07-Mar-2016 Pruning of ISP
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCT0   , ONLY : LSMSSIG, LECMWF
USE YOMMP0   , ONLY : MYPROC
USE YOMRIP   , ONLY : TRIP
USE YOMRES   , ONLY : JPNWST, NRESTS, N1RFS, N2RFS, NFRRES, &
 & NFLSTOP, NFLREST, LSDHM, LDELRES

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: J

REAL(KIND=JPRB) :: ZUNIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "posnam.intfb.h"
#include "namres.nam.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURES',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!        1.1 Set implicit default values

!     ZUNIT  : NUMBER OF SECONDS IN TIME UNIT
ZUNIT=3600._JPRB

NFRRES=-48
N1RFS  =1
N2RFS  =1
DO J=0,JPNWST
  NRESTS(J)=0
ENDDO
LSDHM=.TRUE.
LDELRES=.FALSE.

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
  LSDHM=.FALSE.
  LDELRES=.TRUE.
  NFRRES=1
  NRESTS(0)=-11
  NRESTS(1)=-6
  NRESTS(2)=-12
  NRESTS(3)=-18
  NRESTS(4)=-24
  NRESTS(5)=-30
  NRESTS(6)=-36
  NRESTS(7)=-42
  NRESTS(8)=-48
  NRESTS(9)=-54
  NRESTS(10)=-60
  NRESTS(11)=-66
ENDIF
!      ----------------------------------------------------------------

!*       2.    Modifies default values.
!              ------------------------

CALL POSNAM(NULNAM,'NAMRES')
READ(NULNAM,NAMRES)

!  allow for restart frequency to be in hours
IF(NFRRES < 0.AND.TSTEP > 0.0_JPRB) THEN
  NFRRES= (REAL(-NFRRES,JPRB)*ZUNIT+0.5_JPRB)/TSTEP
ENDIF

!      -----------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' NRESTS = '',40(1X,I0,1X))')NRESTS
WRITE(UNIT=KULOUT,FMT='('' NFRRES = '',1X,I6)')NFRRES
WRITE(UNIT=KULOUT,FMT='('' LSDHM = '',L2)') LSDHM
WRITE(UNIT=KULOUT,FMT='('' LDELRES = '',L2)') LDELRES

!      -----------------------------------------------------------

!*       4.    Call signal handler 
!              -------------------

NFLSTOP=0
NFLREST=0
IF(LSMSSIG) CALL IFSSIG(NFLSTOP,NFLREST,MYPROC)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURES',1,ZHOOK_HANDLE)
END SUBROUTINE SURES
