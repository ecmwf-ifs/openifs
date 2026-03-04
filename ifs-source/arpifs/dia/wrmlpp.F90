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

SUBROUTINE WRMLPP(YDGEOMETRY,YDGFL,YDGFL5,YDSURF,YDSPEC,YDCFU,YDXFU,YDMODEL,CDCONF,PTRAJEC,YDMCUF,YDECV)

!**** *WRMLPP*  - writes out the model level fields of the model

!     Purpose.
!     --------
!     Write out the model level fields of the model

!**   Interface.
!     ----------
!        *CALL* *WRMLPP(...)

!        Explicit arguments :     CDCONF - configuration of call
!        --------------------

!        Implicit arguments :      None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!        Note de travail ARPEGE NR 17

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-03-02

!     Modifications.
!     --------------
!      R. El Khatib 02-03-19 lagged Fourier diffusion before writing out
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad (Dec 2003): cleaning of horizontal diffusion.
!      M.Hamrud      01-DEC-2003 CY28R1 Cleaning
!      R. El Khatib : 01-Mar-2012 LFPOS => LECFPOS
!      F. Vana  28-Nov-2013 : Redesigned trajectory handling
!      R. El Khatib : 23-Aug-2016 lgrbop <=> larpegef
!      P. Lopez: Dec 2015  Write output fields in TL evolution test even if LECFPOS=T
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGFL             , ONLY : TGFL
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : LARPEGEF, LECFPOS
USE YOMTRAJ            , ONLY : TRAJ_TYPE
USE YOMMCUF            , ONLY : TMCUF
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE TYPE_ECV           , ONLY : ECV_CONTAINER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL5
TYPE(TSURF)                  , INTENT(INOUT) :: YDSURF
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSPEC
TYPE(TCFU)                   , INTENT(INOUT) :: YDCFU
TYPE(TXFU)                   , INTENT(INOUT) :: YDXFU
TYPE(MODEL)                  , INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=1)             , INTENT(IN)    :: CDCONF 
TYPE(TRAJ_TYPE)    , OPTIONAL, INTENT(IN)    :: PTRAJEC
TYPE(TMCUF)        , OPTIONAL, INTENT(INOUT) :: YDMCUF
TYPE(ECV_CONTAINER), OPTIONAL, INTENT(IN)    :: YDECV

!     ------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "wrmlppa.intfb.h"
#include "wrmlppg.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRMLPP',0,ZHOOK_HANDLE)
ASSOCIATE(YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI)
ASSOCIATE(LTLEVOL=>YDEPHLI%LTLEVOL)
!     ------------------------------------------------------------------

!*       1.    CALL APPROPRIATE OUTPUT ROUTINE
!              -------------------------------

IF (LARPEGEF) THEN

!*       1.1   ARPEGE FILE OUTPUT

  CALL WRMLPPA(YDGEOMETRY,YDGFL,YDSURF,YDSPEC,YDCFU,YDXFU,YDMODEL,CDCONF,YDMCUF=YDMCUF)

ELSEIF (.NOT.LECFPOS.OR.(LECFPOS.AND.LTLEVOL)) THEN

!*       1.2   GRIB OUTPUT

  CALL WRMLPPG(YDGEOMETRY,YDGFL,YDGFL5,YDSURF,YDMODEL,YDSPEC,PTRAJEC=PTRAJEC,YDECV=YDECV)

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRMLPP',1,ZHOOK_HANDLE)
END SUBROUTINE WRMLPP
