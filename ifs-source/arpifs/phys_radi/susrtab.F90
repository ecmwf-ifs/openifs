! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUSRTAB

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** AER'S RRTM SW RADIATION **

!     M.J. IACONO, AER 
!     JJMorcrette, ECMWF, 20080724: adapted to IFS

!     -----------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESRTAB , ONLY : TRANS, BPADE, RODLOW, RTBLINT

IMPLICIT NONE

INTEGER(KIND=JPIM) :: ITR

REAL(KIND=JPRB) :: ZTAU, ZTFN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUSRTAB',0,ZHOOK_HANDLE)

RODLOW = 0.06_JPRB
RTBLINT= 10000._JPRB
BPADE = 1.0_JPRB / 0.278_JPRB
TRANS(0)   =1.0_JPRB
TRANS(10000)=0.0_JPRB
DO ITR=1,9999
  ZTFN=REAL(ITR)/10000._JPRB
  ZTAU=BPADE*ZTFN/(1.0_JPRB-ZTFN)
  TRANS(ITR)=EXP(-ZTAU)
ENDDO

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSRTAB',1,ZHOOK_HANDLE)
END SUBROUTINE SUSRTAB
