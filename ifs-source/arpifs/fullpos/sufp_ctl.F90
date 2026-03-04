! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFP_CTL

!**** *SUFP_CTL*  - INITIALIZE FULL POST PROCESSING MODULE

!     PURPOSE.
!     --------
!        Initialize variables of post-porcessing shared by all post-processors objects

!**   INTERFACE.
!     ----------
!       *CALL* *SUFP_CTL(...)

!        EXPLICIT ARGUMENTS  
!         -------------------

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 06-Nov-2017 - reconversion of the former subroutine

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "sufpc.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFP_CTL',0,ZHOOK_HANDLE)

CALL SUFPC(LDMODULE=.TRUE.,LDPRINT=.TRUE.)

IF (LHOOK) CALL DR_HOOK('SUFP_CTL',1,ZHOOK_HANDLE)
END SUBROUTINE SUFP_CTL
