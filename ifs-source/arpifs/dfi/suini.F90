! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUINI

!**** *SUINI*  - INITIALIZE MODULE YOMINI

!     PURPOSE
!     -------
!      INITIALIZE CONTROL VARIABLES FOR INITIALIZATION

!**   INTERFACE
!     ---------
!        *CALL* *SUINI*

!     EXTERNALS: None
!     ----------

!     EXPLICIT ARGUMENTS : None
!     ------------------

!     REFERENCE
!     ---------
!     ARPEGE/ALADIN documentation

!     AUTHOR
!     ------
!      Gabor Radnoti GMAP/MICECO
!      Original: 1992.12.24

!     MODIFICATIONS
!     -------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P.Termonia   : 08-12-22 move the call to SUDFI to SU0YOMB (for SSDFI) before cy37
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!     -------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMINI   , ONLY : LDFI, LBIAS, LINCR, LINITER, LSCRINI

!     -------------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

!#include "abor1.intfb.h"
#include "namini.nam.h"
#include "posnam.intfb.h"

!     -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUINI',0,ZHOOK_HANDLE)
!     -------------------------------------------------------------------------

!*      1. SET UP DEFAULT VALUES FOR YOMINI
!          --------------------------------

!        1.1 Set implicit default values

LDFI=.FALSE.
LBIAS=.FALSE.
LINCR=.FALSE.
LINITER=.FALSE.
LSCRINI=.FALSE.

!     -------------------------------------------------------------------------

!*      2. READ NAMINI
!          -----------

CALL POSNAM(NULNAM,'NAMINI')
READ(NULNAM,NAMINI)

IF (.NOT.LDFI) THEN
  LBIAS= .FALSE.
  LINCR= .FALSE.
ENDIF
IF (LBIAS) THEN
  LINCR= .FALSE.
ENDIF

!     -------------------------------------------------------------------------

!*      3. PRINT FINAL VALUES AND WARNING
!          -------------------------------

WRITE(NULOUT,*) '****Module /YOMINI/****'
WRITE(NULOUT,*) 'LDFI = ', LDFI

IF (LDFI) THEN
  WRITE(NULOUT,*) 'INITIAL FIELDS ARE INITIALIZED BY DFI IN ARPEGE/ALADIN'
ELSE
  WRITE(NULOUT,*) 'INITIAL FIELDS ARE NOT INITIALIZED IN ARPEGE/ALADIN'
ENDIF
IF (LBIAS) THEN
  WRITE(NULOUT,*) 'COMPUTING INITIALIZATION INCREMENT ; NSTOP RESET TO 0 !'
ENDIF
IF (LINCR) WRITE(NULOUT,*) 'INCREMENTAL INITIALIZATION'

!     -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUINI',1,ZHOOK_HANDLE)
END SUBROUTINE SUINI
