MODULE SUSTHF_MOD
CONTAINS
SUBROUTINE SUSTHF(YDCST)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_THF  , ONLY : R4LES, R5LES, RHOH2O, R2ES, R4IES, R3LES, R3IES, R5IES, RVTMP2

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *SUSTHF* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOS_CST*
!              This contains the fundamental model constants

!     INTERFACE.
!     ----------
!     CALL *SUSTHF* FROM *SUSURF*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     Original    P. Viterbo      May 2005
!        Based on the IFS suphec


IMPLICIT NONE

TYPE(TCST), INTENT(IN) :: YDCST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSTHF_MOD:SUSTHF',0,ZHOOK_HANDLE)
ASSOCIATE(RATM=>YDCST%RATM, RD=>YDCST%RD, RTT=>YDCST%RTT, RV=>YDCST%RV)

!  DEFINE FUNDAMENTAL CONSTANTS.
!RVTMP2=RCPV/RCPD-1.0_JPRB   !use cp,moist
RVTMP2=0.0_JPRB              !neglect cp,moist
RHOH2O=RATM/100._JPRB
R2ES=611.21_JPRB*RD/RV
R3LES=17.502_JPRB
R3IES=22.587_JPRB
R4LES=32.19_JPRB
R4IES=-0.7_JPRB
R5LES=R3LES*(RTT-R4LES)
R5IES=R3IES*(RTT-R4IES)


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSTHF_MOD:SUSTHF',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUSTHF
END MODULE SUSTHF_MOD
