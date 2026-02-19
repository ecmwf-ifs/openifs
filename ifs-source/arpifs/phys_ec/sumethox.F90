! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUMETHOX

!**** *SUCLDP*   - INITIALIZE MODULE YOEMETH CONTROLLING *METHOX*

!     PURPOSE.
!     --------
!           INITIALIZE YOEMETH

!**   INTERFACE.
!     ----------
!        CALL *SUMETHOX* FROM *SUPHEC*
!              --------        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        MODULE YOEMETH

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

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-04-07
!        Modified : 02-01-29  A.Simmons: increase RQLIM from 3.75e-6
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R.Forbes   01-Mar-2017 increased RQLIM from 4.25e-6 to 4.81e-6
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEMETH   , ONLY : RALPHA1 ,RALPHA2  ,RQLIM   ,&
 & RPBOTOX,  RPBOTPH ,RPTOPOX  ,RPTOPPH ,&
 & RALPHA3,  RLOGPPH  
USE YOMLUN   , ONLY : NULOUT, NULNAM

IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "nammethox.nam.h"
#include "posnam.intfb.h"

!*       1.    SET VALUES
!              ----------

IF (LHOOK) CALL DR_HOOK('SUMETHOX',0,ZHOOK_HANDLE)
RALPHA1=(19._JPRB*LOG(10._JPRB))/(LOG(20._JPRB)**4)
RALPHA2=LOG(1.0_JPRB/3._JPRB+0.01_JPRB)
RQLIM  =4.81E-6_JPRB    ! relaxation limit equivalent to 7.7ppmv
RPBOTOX=10000._JPRB
RPBOTPH=20._JPRB
RPTOPOX=50._JPRB
RPTOPPH=0.1_JPRB
RALPHA3=0.5_JPRB*(LOG(100._JPRB)+RALPHA2)
RLOGPPH=LOG(RPTOPPH/RPBOTPH)

CALL POSNAM(NULNAM,'NAMMETHOX')
READ(NULNAM,NAMMETHOX)
WRITE(UNIT=NULOUT,FMT='('' COMMON YOEMETH '')')
WRITE(UNIT=NULOUT,FMT='('' RQLIM= '',E12.5)') RQLIM

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUMETHOX',1,ZHOOK_HANDLE)
END SUBROUTINE SUMETHOX
