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

SUBROUTINE SUALPGO

!**** *SUALPGO*  - DDH output : allocate space for pseudo-grib coding

!     Purpose.
!     --------
!     Allocate space for real and integer blocks for subsequent
!     pseudo grib coding and output of DDH data.

!**   Interface.
!     ----------
!        *CALL* *SUALPGO(..)*

!        Explicit arguments :     None
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

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-03-01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LALLOPR
USE YOMPGO   , ONLY : NLENREA  ,NLENIND  ,NLENINT
USE YOMPGOM  , ONLY : REAPG    ,NINDPG   ,NINTPG
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.  ALLOCATE SPACE.
!            ---------------

IF (LHOOK) CALL DR_HOOK('SUALPGO',0,ZHOOK_HANDLE)
ALLOCATE(REAPG (NLENREA))
ALLOCATE(NINDPG(NLENIND))
ALLOCATE(NINTPG(NLENINT))
IF (LALLOPR) THEN
  WRITE(NULOUT,9990) 'REAPG     ',SIZE(REAPG),SHAPE(REAPG)
  WRITE(NULOUT,9990) 'NINDPG    ',SIZE(NINDPG),SHAPE(NINDPG)
  WRITE(NULOUT,9990) 'NINTPG    ',SIZE(NINTPG),SHAPE(NINTPG)
ENDIF

IF (LHOOK) CALL DR_HOOK('SUALPGO',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9990 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE SUALPGO
