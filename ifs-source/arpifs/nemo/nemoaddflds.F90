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

SUBROUTINE NEMOADDFLDS(YDMCC,KSTGLO,KIDIA,KFDIA,PSST,PTSK)
!
!**** *NEMOADDFLDS*  - Accumulate flux data to be sent to the ocean and sea ice models
!
!     Purpose.
!     --------
!       UPDATE DIAGNOSTICS FIELDS TO BE SENT TO THE OCEAN
!
!**   Interface.
!     ----------
!       *CALL*  *NEMOADDFLDS(KSTGLO, KIDIA, KFDIA, PSST, PTSK)
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
!       
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!
!     Modifications.
!     --------------
!
!     -----------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMCC   , ONLY : TMCC
USE CPLNG

IMPLICIT NONE

TYPE(TMCC)         ,INTENT(INOUT):: YDMCC
INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO, KIDIA, KFDIA
REAL(KIND=JPRB),    INTENT(IN), DIMENSION(KIDIA:KFDIA) :: PSST, PTSK

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IST,IEND

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NEMOADDFLDS',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & LNEMOCOUP=>YDMCC%LNEMOCOUP, CPLNG_FLD=>YDMCC%CPLNG_FLD &
 & )

IF (LNEMOCOUP) THEN

   IST  = KSTGLO - 1 + KIDIA
   IEND = KSTGLO - 1 + KFDIA
   CPLNG_FLD(YDMCC%IP_A_SST_ATM)%D(IST:IEND,1,1) = PSST(:)
   CPLNG_FLD(YDMCC%IP_A_TSK_ATM)%D(IST:IEND,1,1) = PTSK(:)

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NEMOADDFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE NEMOADDFLDS
