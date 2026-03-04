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

SUBROUTINE SUSPE0(YDGEOMETRY,YDDIMF,YDSPEC)

!**** *SUSPE0*  - Initialize the spectral fields to 0.

!     Purpose.
!     --------
!           Initialize spectral fields to 0.

!**   Interface.
!     ----------
!        *CALL* *SUSPE0*

!        Explicit arguments :
!        --------------------
!        YDGEOMETRY
!        YDSPEC  (replaces the previously implicit argument SPA3)


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
!      Original : 87-12-14

!     Modifications.
!     --------------
!      O.Marsden     August 2016 Work on explicit spectral field argument
!      K. Yessad (Dec 2016): Prune obsolete options.
!      -------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDIMF  , ONLY : TDIMF
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!      -----------------------------------------------------------------
IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT(IN)    :: YDGEOMETRY
TYPE(TDIMF)          ,INTENT(INOUT) :: YDDIMF
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSPEC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSPE0',0,ZHOOK_HANDLE)
ASSOCIATE(YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NFD2D=>YDDIMF%NFD2D,NPSP=>YDMP%NPSP)
!      -----------------------------------------------------------------

!*       1.    Initialize spectral fields
!              ---------------------------

YDSPEC%SP3D(:,:,:)=0.0_JPRB

IF(NPSP == 1)THEN
  YDSPEC%SP2D(:,1:NFD2D)=0.0_JPRB
ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSPE0',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPE0
