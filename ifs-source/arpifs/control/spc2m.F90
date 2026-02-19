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

SUBROUTINE SPC2M(YDGEOMETRY,YDRIP,YDDYN,CDCONF,YDSPEC)

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE YOMRIP             , ONLY : TRIP
USE YOMDYN             , ONLY : TDYN
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!**** *SPC2M* - Multi tasking interface to SPC2

!     Purpose.
!     --------
!        Multi tasking interface (also calls SPC2 single-tasked if
!        multi-tasking not requested).

!**   Interface.
!     ----------
!        *CALL* *SPC2M(...)

!        Explicit arguments :  CDCONF - configuration of work (see doc.)

!        Implicit arguments :
!        --------------------
!                            NONE.
!     Method.
!     -------

!     Externals.   SPC2 - Spectral space computations - 2-D case
!     ----------   Called by SPCH.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-11-24

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      O. Marsden       May 2016  Remove redundant geometry argument
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)          ,INTENT(INOUT) :: YDDYN
TYPE(TRIP)          ,INTENT(INOUT) :: YDRIP
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF 
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JM, IM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "spc2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPC2M',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NUMP=>YDDIM%NUMP)
!     ------------------------------------------------------------------

!*       1.    LOOP OVER ZONAL WAVENUMBER M.
!              -----------------------------

DO JM=1,NUMP
  IM=YDLAP%MYMS(JM)
  CALL SPC2(YDGEOMETRY,YDRIP,YDDYN,CDCONF,IM,JM,YDSPEC)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPC2M',1,ZHOOK_HANDLE)
END SUBROUTINE SPC2M
