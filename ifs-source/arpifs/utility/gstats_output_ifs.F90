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

SUBROUTINE GSTATS_OUTPUT_IFS(YDRIP)

!**** *GSTATS_OUTPUT_IFS* - print timing statistics

!     PURPOSE.
!     --------
!       To print out timings gathered by GSTATS

!**   INTERFACE.
!     ----------
!       *CALL* *GSTATS_OUTPUT_IFS*

!        EXPLICIT ARGUMENTS     None
!        --------------------

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.   
!     ----------   

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      Mats Hamrud ECMWF
!      ORIGINAL : 98-11-15

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      M.Drusch      18-Jan-2007 introduce nconf 302
!      R. El Khatib  26-Feb-2008 Bugfix
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPROC
USE YOMCT0   , ONLY : NCONF
USE YOMRIP   , ONLY : TRIP
USE YOMGSTATS, ONLY : LSTATS

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP),INTENT(INOUT):: YDRIP
REAL(KIND=JPRD) :: ZAVEAVE(0:100)
REAL(KIND=JPRB) :: ZMEAN,ZFDPD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GSTATS_OUTPUT_IFS',0,ZHOOK_HANDLE)
ASSOCIATE(NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------

IF (LSTATS) THEN
  CALL GSTATS_PRINT(NULOUT,ZAVEAVE,100)
  ZMEAN = REAL(ZAVEAVE(1),JPRB)/(1000._JPRB*NPROC)
  IF( (NCONF == 1.OR.NCONF == 302) .AND. ZMEAN > 0.0_JPRB )THEN
    ZFDPD=REAL(NSTOP,JPRB)*TSTEP/ZMEAN
    WRITE(NULOUT,'(A,F10.1)')'FORECAST DAYS PER DAY ',ZFDPD
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GSTATS_OUTPUT_IFS',1,ZHOOK_HANDLE)
END SUBROUTINE GSTATS_OUTPUT_IFS
