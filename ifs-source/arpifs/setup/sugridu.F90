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

SUBROUTINE SUGRIDU(YDGEOMETRY,YDEPHLI,YDML_GCONF,YDML_LBC,YDSPEC,YDGFL,KFILE)

!**** *SUGRIDU*  - Routine to initialize the upper air grid point
!                  fields of the model.

!     Purpose.
!     --------
!           Initialize the upper air grid point fields of the model.

!**   Interface.
!     ----------
!        *CALL* *SUGRIDU*

!        Explicit arguments :
!        --------------------
!        KFILE : an indicator for which  file is to be read

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  SUGRIDUG - initialize from GRIB file
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 94-06-30

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!     ------------------------------------------------------------------

USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOEPHLI                , ONLY : TEPHLI
USE YEMLBC_MODEL             , ONLY : TELBC_MODEL
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPIM     ,JPRB
USE YOMHOOK                , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCT0                 , ONLY : N3DINI
USE YOMARG                 , ONLY : NGRIBFILE
USE SPECTRAL_FIELDS_MOD    , ONLY : ASSIGNMENT(=), SPECTRAL_FIELD

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TEPHLI)                 , INTENT(INOUT) :: YDEPHLI
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(INOUT) :: YDML_GCONF
TYPE(TELBC_MODEL)                , INTENT(IN)    :: YDML_LBC
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSPEC
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KFILE 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "sugridua.intfb.h"
#include "sugridug.intfb.h"
#include "sugridug1.intfb.h"
#include "sugridug2.intfb.h"

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGRIDU',0,ZHOOK_HANDLE)
!      -----------------------------------------------------------

CALL GSTATS(24,0)

!*       1.    INITIALIZE UPPER AIR GRIDPOINT FIELD.
!              -------------------------------------

IF (NGRIBFILE /= 1) THEN

!*      1.1   INITIALIZE FIELDS FROM *FA* FILE.

  CALL SUGRIDUA(YDGEOMETRY,YDGFL,YDML_GCONF,YDML_LBC,KFILE)

ELSE

!*      1.2   INITIALIZE FIELDS FROM GRIB FILES OR ARTIFICIAL DATA

  IF ( N3DINI == 0 ) THEN
    CALL SUGRIDUG(YDGEOMETRY,YDML_GCONF,YDGFL,YDEPHLI,YDML_GCONF%YGFL,KFILE)
  ELSEIF( N3DINI == 1 ) THEN
    CALL SUGRIDUG1(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL)
  ELSEIF( N3DINI == 2 ) THEN
    CALL SUGRIDUG2(YDGEOMETRY,YDML_GCONF,YDSPEC,YDGFL)
  ELSE
    CALL ABOR1('SUGRIGU: ONLY N3DINI = 0, 1 or 2 are implemented currently')
  ENDIF

ENDIF

!     ------------------------------------------------------------------

CALL GSTATS(24,1)

IF (LHOOK) CALL DR_HOOK('SUGRIDU',1,ZHOOK_HANDLE)
END SUBROUTINE SUGRIDU
