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

SUBROUTINE SUIOS

!**** *SUIOS* - ROUTINE TO SET UP FOR USING WORK FILES

!     Purpose.   TO SET UP COMMMON BLOCK YOMIOS WHICH CONTAINS   -
!     --------   PARAMETERS FOR USING THE MIO PACKAGE ON WORK
!                FILES. OPENS WORK FILES.

!**   Interface.
!     ----------
!        *CALL* *SUIOS*

!        Explicit arguments :  None.
!        --------------------

!        Implicit arguments : COMMON YOMIOS
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud   *ECMWF*
!      Original : 91-01-28

!     Modifications.
!     --------------
!      Modified : 01-06-13 R. El Khatib Restart file name for ISP
!      Modified : 02-06-10 D.Dent  Restart file name for fields
!                                  which are not part of the surface workfile
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad : 04-07-06 Old SUEIOS now in SUIOS.
!      K. Yessad (Jan 2010): remove useless variables.
!      R. El Khatib 07-Mar-2016 Pruning of ISP
!     ------------------------------------------------------------------


USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : CNMEXP, LECMWF
USE YOMIOS   , ONLY : CIOSPRF, CFRCF, LRCFTIME
USE YOMLUN   , ONLY : NULNAM

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namios.nam.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUIOS',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    PARAMETERS UNDER USER CONTROL.
!              ------------------------------

!*       1.1   SET DEFAULT VALUES

!        1.1.1 Set implicit default values

CFRCF    = 'rcf'
CIOSPRF  = 'srf'
LRCFTIME = .FALSE.

!        1.1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
  CFRCF   = 'RFILE'//CNMEXP(1:4)
  CIOSPRF = 'RSPEC'//CNMEXP(1:4)
ENDIF

!*       1.2   READ NAMELIST.

CALL POSNAM(NULNAM,'NAMIOS')
READ(NULNAM,NAMIOS)

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUIOS',1,ZHOOK_HANDLE)
END SUBROUTINE SUIOS
