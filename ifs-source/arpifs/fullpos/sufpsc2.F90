! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPSC2(KFPROMA,KFPWIDE)

!**** *SUFPSC2*  - FULL-POS HORIZONTAL SCANNING

!     PURPOSE.
!     --------
!         INITIALIZE CACHE-BLOCKING FACTOR FOR POST-PROCESSING ARRAYS IN THE TARGET GEOMETRY
!         INITIALIZE WIDTH OF HALOS FOR INTERPOLATIONS (not necessary, actually)

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPSC2*

!        EXPLICIT ARGUMENTS
!        ------------------
!        KFPROMA : default segmentation
!        KFPWIDE : width of halos

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
!      RYAD EL KHATIB *METEO-FRANCE*
!      G. Radnoti   : 95-01-23   full pos of aladin

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-01-29 Remove SM aspects/cleanings
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 23-may-2005 NFPBLOFF
!      M. Jidane    : 19-04-2006  Correction of a bug in allocation of NFPBLOFF
!      K. Yessad    : 27-Feb-2007 Adapt to arrival geometry, clean.
!      R. El Khatib : 16-Jul-2012 Enable NFPROMA < 0
!      R. El Khatib 27-Jul-2016 setup of NFPSLWIDE + code reordering
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULNAM

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(OUT) :: KFPWIDE
INTEGER(KIND=JPIM), INTENT(INOUT) :: KFPROMA

INTEGER(KIND=JPIM) :: NFPROMA, NFPSLWIDE ! Sorry Doctor : namelist variables

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namfpsc2.nam.h"

!     ------------------------------------------------------------------

#include "posnam.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPSC2',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1. COMPUTES SEGMENTATION 
!           ---------------------

NFPROMA=KFPROMA
NFPSLWIDE=0

!      1.1 Read namelist & control

CALL POSNAM(NULNAM,'NAMFPSC2')
READ(NULNAM,NAMFPSC2)

KFPWIDE=NFPSLWIDE
KFPROMA=NFPROMA

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPSC2',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPSC2
