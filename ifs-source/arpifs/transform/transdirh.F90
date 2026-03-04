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

SUBROUTINE TRANSDIRH(YDGEOMETRY,LDGRADSP,YDGFL,YDGMV,CDCONF,KNFTHER,YDSP)

!**** *TRANSDIRH * - Direct transforms

!     Purpose.  Perform direct transform (gridpoint to spectral)
!     --------

!     Explicit arguments :
!     --------------------
!        CDCONF     - configuration of work
!                     possible values: A, G or U

!     Externals.
!     ----------
!     TRANSDIR_MDL - transforms for thr model   

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-10-25
!        Modified : 03-08-01 M.Hamrud - GFL introduction
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!        R. El Khatib 18-Jul-2012 Fullpos move away
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL       , ONLY : TGFL
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPRB
USE YOMHOOK      , ONLY : LHOOK,   DR_HOOK, JPHOOK
!USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE SPECTRAL_FIELDS_MOD

IMPLICIT NONE
TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
LOGICAL             ,INTENT(IN)    :: LDGRADSP
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNFTHER
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "transdir_mdl.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSDIRH',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CALL TRANSDIR_MDL(YDGEOMETRY,LDGRADSP,YDGFL,YDGMV,CDCONF,KNFTHER,YDSP)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSDIRH',1,ZHOOK_HANDLE)
END SUBROUTINE TRANSDIRH

