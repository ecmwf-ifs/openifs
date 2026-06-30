SUBROUTINE SUVDF
! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *SUVDF* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEVDF*

!     A.C.M. BELJAARS         E.C.M.W.F.       2/11/89

!     PURPOSE
!     -------

!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOEVDF*

!     INTERFACE.
!     ----------

!     CALL *SUVDF* FROM *SUPHEC*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     NONE.

!     REFERENCE.
!     ----------

!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB         ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEVDF   , ONLY : RLAM     ,RKAP     ,RVDIFTS  ,&
 & REPDU2   ,RENTR    ,&
 & RPAR     ,RPAR1    ,RPARSRF 

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.     SET FIRST SET OF CONSTANTS
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('SUVDF',0,ZHOOK_HANDLE)

RLAM   =150._JPRB
RKAP   =0.4_JPRB
RVDIFTS=1.0_JPRB

!     ------------------------------------------------------------------

!*         2.      SET OTHER CONSTANTS
!                  -------------------

REPDU2 =(0.1_JPRB)**2

!     ENTRAINMENT PARAMETRIZATION

RENTR=0.20_JPRB
RPAR=2._JPRB
RPAR1=0.6_JPRB
RPARSRF=0.1_JPRB

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUVDF',1,ZHOOK_HANDLE)
END SUBROUTINE SUVDF
