! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_MCICA(YDEMCICA)

!**** *SU_MCICA*   - INITIALIZE COMMON YOE_MCICA

!     PURPOSE.
!     --------
!           INITIALIZE YOE_MCICA, THE COMMON THAT CONTAINS COEFFICIENTS
!           NEEDED TO RUN THE RAISANEN-COLE-BARKER VERSION OF McICA

!**   INTERFACE.
!     ----------
!        *CALL* *SU_MCICA* FROM *SUECRAD*


!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOE_MCICA


!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE *ECMWF*
!        from Jason Cole's XCW_DATA.dat

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2006-01-05
!     G.Mozdzynski March 2011 read constants from files
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB, JPIM
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOE_MCICA, ONLY : TEMCICA
USE YOMLUN    ,ONLY : RESERVE_LUN, FREE_LUN
USE YOMMP0    ,ONLY : NPROC, MYPROC
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD

IMPLICIT NONE

TYPE(TEMCICA),INTENT(INOUT):: YDEMCICA
INTEGER(KIND=JPIM) :: IULTMP

CHARACTER(LEN = 512) :: CLZZZ
CHARACTER(LEN = 512) :: CLF1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU_MCICA',0,ZHOOK_HANDLE)
ASSOCIATE(NMCI1=>YDEMCICA%NMCI1, NMCI2=>YDEMCICA%NMCI2, XCW=>YDEMCICA%XCW, XCW_D=>YDEMCICA%XCW_D)
NMCI1 = 1000
NMCI2 = 140

IF( MYPROC==1 )THEN
  IULTMP = RESERVE_LUN()
  CALL GET_ENVIRONMENT_VARIABLE("DATA",CLZZZ)
  IF(CLZZZ /= " ") THEN
    CLF1=TRIM(CLZZZ)//"/ifsdata/MCICA"
    WRITE(0,'(1x,a)')'Reading '//TRIM(CLF1)
#ifdef LITTLE_ENDIAN
    OPEN(IULTMP,FILE=CLF1,FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')
#else
    OPEN(IULTMP,FILE=CLF1,FORM="UNFORMATTED",ACTION="READ",ERR=1000)
#endif
  ELSE
#ifdef LITTLE_ENDIAN
    OPEN(IULTMP,FILE='MCICA',FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')
#else
    OPEN(IULTMP,FILE='MCICA',FORM="UNFORMATTED",ACTION="READ",ERR=1000)
#endif
  ENDIF

  READ(IULTMP,ERR=1001) XCW_D
  CLOSE(IULTMP)
  CALL FREE_LUN(IULTMP)
  XCW = REAL(XCW_D,JPRB)

ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (XCW,MTAGRAD,1,CDSTRING='SU_MCICA:')
ENDIF
  
!     ----------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SU_MCICA',1,ZHOOK_HANDLE)
RETURN

1000 CONTINUE
CALL ABOR1("SU_MCICA:ERROR OPENING FILE MCICA")
1001 CONTINUE
CALL ABOR1("SU_MCICA:ERROR READING FILE MCICA")

END SUBROUTINE SU_MCICA
