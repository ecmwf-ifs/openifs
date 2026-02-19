! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_CLOP550 

!**** *SU_CLOP550*  - INITIALIZE COMMON YOECLOP550

!     PURPOSE.
!     --------
!           INITIALIZE YOECLOP550, WITH CLOUD OPTICAL PARAMETERS

!**   INTERFACE.
!     ----------
!        *CALL*  SUCLOP550

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLOP550

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     A. Slingo, 1989: J. Atmos. Sci., 46, 1419-1427
!     Fu, 1996: J. Climate, 9, 2058-2082

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 20060324
!        Vincent Huijnen : 201301112 - Modification of RSA55

!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOECLOP550, ONLY : RSA55, RSB55, RSC55, RSD55, RFA55

IMPLICIT NONE

!     -----------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_CLOP550',0,ZHOOK_HANDLE)

!*          1.    SHORTWAVE CLOUD OPTICAL PROPERTIES AT 550 NM
!                 --------------------------------------------

!-- WATER CLOUDS FROM SLINGO (1989)

!! try 330 nm values
!RSA55=3.308E-02_JPRB 
!RSB55=1.246_JPRB
!RSC55=-3.0E-7_JPRB
!RSD55=2.36E-7_JPRB
!RFA55 = (/-0.2937E-03_JPRB, 0.2545E+01_JPRB /)


! try 550 nm values
RSA55=2.838E-02_JPRB 
RSB55=1.300_JPRB
RSC55=0._JPRB
RSD55=0._JPRB
!-- ICE CLOUDS FROM FU (1996)

RFA55 = (/-0.303108E-04_JPRB, 0.251805E+01_JPRB /)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU_CLOP550',1,ZHOOK_HANDLE)
END SUBROUTINE SU_CLOP550
