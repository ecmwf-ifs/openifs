! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUECOZC(YDEOZOC)

!**** *SUECOZC* - DEFINES CLIMATOLOGICAL DISTRIBUTION OF OZONE
!                 (FORTUIN LANGEMATZ CLIMATOLOGY OF OZONE)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECOZC* FROM ******** a setup routine somewhere! suecrad?

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===
! ROZT    :                : AMOUNT OF OZONE (KG/KG) 

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     J.-J. MORCRETTE  E.C.M.W.F.    95/01/25

!     MODIFICATIONS.
!     --------------
!     JJMorcrette 971117 Cy18
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     F. Vana  05-Mar-2015  Support for single precision
!     O. Marsden Jan  2017  Moved non-setup code to a new UPDECOZC routine
!-----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPRD, JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEOZOC   , ONLY : TEOZOC
USE YOMLUN    , ONLY : RESERVE_LUN, FREE_LUN
USE YOMMP0    , ONLY : NPROC, MYPROC
USE MPL_MODULE, ONLY : MPL_BROADCAST
USE YOMTAG    , ONLY : MTAGRAD

IMPLICIT NONE

!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------
TYPE(TEOZOC)      ,INTENT(INOUT) :: YDEOZOC

!*       0.2   LOCAL ARRAYS.
!              -------------

REAL(KIND=JPRD),ALLOCATABLE :: ZOZCL_D(:,:,:)

INTEGER(KIND=JPIM) :: IULTMP
CHARACTER(LEN = 256) :: CLZZZ
CHARACTER(LEN = 256) :: CLF1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUECOZC',0,ZHOOK_HANDLE)

!* Ozone climatological mixing ratio (ppmv)
!        from Langematz and Fortuin 

  ALLOCATE(YDEOZOC%ZOZCL(19,34,12))

  IF( MYPROC==1 )THEN
    ALLOCATE(ZOZCL_D(19,34,12))
    IULTMP = RESERVE_LUN()
    CALL GET_ENVIRONMENT_VARIABLE("DATA",CLZZZ)
    IF(CLZZZ /= " ") THEN
      CLF1=TRIM(CLZZZ)//"/ifsdata/ECOZC"
      WRITE(0,'(1x,a)')'Reading '//TRIM(CLF1)
#ifdef LITTLE_ENDIAN
      OPEN(IULTMP,FILE=CLF1,FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')
#else
      OPEN(IULTMP,FILE=CLF1,FORM="UNFORMATTED",ACTION="READ",ERR=1000)
#endif
    ELSE
#ifdef LITTLE_ENDIAN
      OPEN(IULTMP,FILE='ECOZC',FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')
#else
      OPEN(IULTMP,FILE='ECOZC',FORM="UNFORMATTED",ACTION="READ",ERR=1000)
#endif
    ENDIF
    READ(IULTMP,ERR=1001) ZOZCL_D
    CLOSE(IULTMP)
    CALL FREE_LUN(IULTMP)
    YDEOZOC%ZOZCL = REAL(ZOZCL_D,JPRB)
    DEALLOCATE(ZOZCL_D)
  ENDIF
  IF( NPROC>1 )THEN
    CALL MPL_BROADCAST (YDEOZOC%ZOZCL,MTAGRAD,1,CDSTRING='SUECOZC:')
  ENDIF

 

! PRESSURE LEVELS FOR O3 CLIMATOLOGY

YDEOZOC%RPROC = (/&
 & 0.0_JPRB&
 & ,                    0.3_JPRB,    0.5_JPRB,    0.7_JPRB,      1.0_JPRB&
 & ,     1.5_JPRB,    2.0_JPRB,     3._JPRB,     5._JPRB,     7._JPRB,     10._JPRB&
 & ,    15._JPRB,    20._JPRB,    30._JPRB,    50._JPRB,    70._JPRB,    100._JPRB&
 & ,   150._JPRB,   200._JPRB,   300._JPRB,   500._JPRB,   700._JPRB,   1000._JPRB&
 & ,  1500._JPRB,  2000._JPRB,  3000._JPRB,  5000._JPRB,  7000._JPRB,  10000._JPRB&
 & , 15000._JPRB, 20000._JPRB, 30000._JPRB, 50000._JPRB, 70000._JPRB, 100000._JPRB&
 & , 110000._JPRB&
 & , 0.0_JPRB   ,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,     0.0_JPRB&
 & , 0.0_JPRB   ,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,     0.0_JPRB&
 & , 0.0_JPRB   ,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,     0.0_JPRB&
 & , 0.0_JPRB   ,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,    0.0_JPRB,     0.0_JPRB&
 & , 0.0_JPRB&
 & /)
  

IF (LHOOK) CALL DR_HOOK('SUECOZC',1,ZHOOK_HANDLE)
RETURN

1000 CONTINUE
CALL ABOR1("SUECOZC:ERROR OPENING FILE ECOZC")
1001 CONTINUE
CALL ABOR1("SUECOZC:ERROR READING FILE ECOZC")

END SUBROUTINE SUECOZC
