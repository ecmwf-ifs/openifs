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

SUBROUTINE SUMSC(KULOUT)

!**** *SUMSC * - Routine to initialize machine specific constants

!     Purpose.
!     --------
!           Initialize machine (compiler) constants

!**   Interface.
!     ----------
!        *CALL* *SUMSC(KULOUT)

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        COMMON YOMMSC

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
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 96-03-20
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND2  ,ONLY : JPRH
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMMSC   , ONLY : NINTLEN  ,NREALEN  ,NLOGLEN  ,NDBLLEN

IMPLICIT NONE

!     EXTERNAL INTEGER FUNCTIONS
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

REAL(KIND=JPRB) :: ZZ(2)
INTEGER(KIND=JPIM) :: II(2)
REAL(KIND=JPRH) :: Z_DL(2)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
!      1. INITIALIZE LENGHTS (IN BYTES) OF DIFFERENT TYPES

IF (LHOOK) CALL DR_HOOK('SUMSC',0,ZHOOK_HANDLE)
ZZ(:) = 1.0_JPRB
II(:) = 1
Z_DL(:) = 1.0_JPRH

NINTLEN = 0
NREALEN = 0
NLOGLEN = 0
NDBLLEN = 0

NINTLEN = STORAGE_SIZE(II)/8
NREALEN = STORAGE_SIZE(ZZ)/8
NLOGLEN = 8
!NDBLLEN = STORAGE_SIZE(Z_DL)/8
NDBLLEN = 16
!!NDBLLEN = STORAGE_SIZE(Z_DL)/8

WRITE(KULOUT,'(4(1X,A,I3))')&
 & 'NINTLEN=',NINTLEN,'NREALEN=',NREALEN,&
 & 'NLOGLEN=',NLOGLEN,'NDBLLEN=',NDBLLEN  

IF( NINTLEN == 0 .OR. NREALEN == 0 .OR. NLOGLEN == 0 .OR. NDBLLEN == 0 )THEN
  CALL ABOR1('SUMSC: WORD KIND LENGTH PROBLEM')
ENDIF

IF (LHOOK) CALL DR_HOOK('SUMSC',1,ZHOOK_HANDLE)
END SUBROUTINE SUMSC
