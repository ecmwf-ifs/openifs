! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPSC2_DEP(KFPROMA_DEP)

!**** *SUFPSC2_DEP*  - FULL-POS HORIZONTAL SCANNING

!     PURPOSE.
!     --------
!        SETUP CACHE-BLOCKING FACTOR FOR POST-PROCESSING ARRAYS ON THE 'DEPARTURE' GEOMETRY

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPSC2_DEP*

!        EXPLICIT ARGUMENTS
!        ------------------
!        NONE.

!        IMPLICIT ARGUMENTS
!        --------------------
!        See modules below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION ABOUT FULL-POS

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE* (SUFPSC2)
!      K. Yessad (rename into SUFPSC2_DEP)
!      G. Radnoti   : 95-01-23   full pos of aladin

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-01-29 Remove SM aspects/cleanings
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 23-may-2005 NFPBLOFF
!      M. Jidane    : 19-04-2006  Correction of a bug in allocation of NFPBLOFF
!      K. Yessad    : 27-Feb-2007 Rename into SUFPSC2_DEP, adapt, clean.
!      R. El Khatib : 17-Sep-2007 Proper default value for NFPROMA_DEP if the default is negative
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2 (differentiation
!      of interpolation grid and output grid)
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULNAM

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(INOUT) :: KFPROMA_DEP

INTEGER(KIND=JPIM) :: NFPROMA_DEP ! sorry Doctor : namelist variable

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namfpsc2_dep.nam.h"

#include "posnam.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPSC2_DEP',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

! Default value
NFPROMA_DEP=KFPROMA_DEP

!      1.2 Read namelist

CALL POSNAM(NULNAM,'NAMFPSC2_DEP')
READ(NULNAM,NAMFPSC2_DEP)

KFPROMA_DEP=NFPROMA_DEP
! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPSC2_DEP',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPSC2_DEP
