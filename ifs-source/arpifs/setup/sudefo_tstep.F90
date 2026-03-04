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

SUBROUTINE SUDEFO_TSTEP(YDDIM,YDDYNA,YDRIP)

!------------------------------------------------------------------------------
!**** *SUDEFO_TSTEP*   - Initialize default value for timestep

!     Purpose.
!     --------
!           Initialize default value for TSTEP

!**   Interface.
!     ----------
!        *CALL* *SUDEFO_TSTEP

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

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
!      Original : 87-10-15 (SUDYN)

! Modifications
! -------------
!   K. Yessad (Jan 2014): put code in SUDEFO_TSTEP.
!   O. Marsden(Nov 2017): removed UTSTEP from YOMARG
! End Modifications
!------------------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMARG   , ONLY : NSUPERSEDE
USE YOMCT0   , ONLY : NCONF
USE YOMRIP   , ONLY : TRIP
USE YOMDYNA  , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN)   :: YDDIM
TYPE(TDYNA), INTENT(IN)   :: YDDYNA
TYPE(TRIP) , INTENT(INOUT):: YDRIP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDEFO_TSTEP',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------


IF (NSUPERSEDE==0 ) THEN
  IF(NSMAX == 639) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=600._JPRB
      ELSE
        TSTEP=300._JPRB
      ENDIF
    ELSE
      TSTEP=100._JPRB
    ENDIF
  ELSEIF(NSMAX == 319) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=1200._JPRB
      ELSE
        TSTEP=600._JPRB
      ENDIF
    ELSE
      TSTEP=300._JPRB
    ENDIF
  ELSEIF(NSMAX == 213) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=1200._JPRB
      ELSE
        TSTEP=600._JPRB
      ENDIF
    ELSE
      TSTEP=300._JPRB
    ENDIF
  ELSEIF(NSMAX == 159) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=3600._JPRB
      ELSE
        TSTEP=1800._JPRB
      ENDIF
    ELSE
      TSTEP=600._JPRB
    ENDIF
  ELSEIF(NSMAX == 106) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=3600._JPRB
      ELSE
        TSTEP=1800._JPRB
      ENDIF
    ELSE
      TSTEP=900._JPRB
    ENDIF
  ELSEIF(NSMAX == 95) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=3600._JPRB
      ELSE
        TSTEP=1800._JPRB
      ENDIF
    ELSE
      TSTEP=1350._JPRB
    ENDIF
  ELSEIF(NSMAX == 63) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=3600._JPRB
      ELSE
        TSTEP=1800._JPRB
      ENDIF
    ELSE
      TSTEP=1350._JPRB
    ENDIF
  ELSEIF(NSMAX < 63) THEN
    IF(YDDYNA%LSLAG) THEN
      IF (YDDYNA%LTWOTL) THEN
        TSTEP=3600._JPRB
      ELSE
        TSTEP=1800._JPRB
      ENDIF
    ELSE
      TSTEP=1800._JPRB
    ENDIF
  ENDIF
  IF (NCONF > 899) TSTEP=1._JPRB
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDEFO_TSTEP',1,ZHOOK_HANDLE)
END SUBROUTINE SUDEFO_TSTEP
