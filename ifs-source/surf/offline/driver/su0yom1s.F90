SUBROUTINE SU0YOM1S
USE PARKIND1    ,ONLY : JPIM     ,JPRB,    JPRD
USE YOMLUN1S    ,ONLY : NULOUT   ,NULNAM
USE YOMRIP      ,ONLY : NINDAT   ,NSSSSS
USE MPL_MODULE
USE YOMHOOK     ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SU0YOM1S*  - INITIALIZE LEVEL 0 COMMONS

!     PURPOSE.
!     --------
!           INITIALIZE LEVEL 0 COMMONS (CONSTANT ALONG ALL THE JOB).

!**   INTERFACE.
!     ----------
!        *CALL* *SU0YOM1S*

!        EXPLICIT ARGUMENTS
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!      SUCST   - INITIALIZE CONSTANTS
!      SU1S    - INITIALIZE LOGICAL SWITCHES
!      SULUN1S - INITIALIZE LOGICAL UNITS
!      SURIP   - INITIALIZE MODEL TIME
!      SUPHEC  - INITIALIZE PHYSICS

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the 1D-surface model

!     AUTHOR.
!     -------
!        Jean Francois Mahfouf and Pedro Viterbo  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL   : 95-03-01
!        BART VD HURK (KNMI): ADJUSTED FOR MULTI-COLUMN MODE
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER(KIND=JPIM) :: NPROC, MYPROC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "su1s.intfb.h"
#include "surip.intfb.h"
#include "sucst.intfb.h"
#include "su0phy1s.intfb.h"
#include "sudim1s.intfb.h"
#include "suoptsurf.intfb.h"
#include "sulun1s.intfb.h"
#include "suphec.intfb.h"

IF (LHOOK) CALL DR_HOOK('SU0YOM1S',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    INITIALIZE COMMON YOMLOG1S
!              --------------------------
NPROC = MPL_NPROC()
MYPROC = MPL_MYRANK()

NULNAM     =  4
OPEN(NULNAM,FILE='input',STATUS='OLD',ACTION='READ')
NULOUT = 6
IF( NPROC > 1 ) THEN
  IF( MYPROC /= 1 ) THEN
    OPEN(UNIT=NULOUT, FILE='/dev/null')
  ENDIF
ENDIF 

CALL SU1S

!*       3.    INITIALIZE MODEL TIME 
!              --------------------- 

CALL SURIP(NULOUT)


!*       4.    INITIALIZE CONSTANTS.
!              ---------------------

CALL SUCST(NULOUT,NINDAT,NSSSSS,0)

!*       5.    INITIALIZE COMMON YOEPHY 
!              ------------------------ 

CALL SU0PHY1S(NULOUT)


!*       6.    INITIALIZE DIMENSIONS.
!              ----------------------

CALL SUDIM1S

!*       7.    INITIALIZE OPTIMIZED SFC PARAMETERS 
!              ------------------------ 

CALL SUOPTSURF(NULOUT)



!*       2.    INITIALIZE COMMON YOMLUN1S (MOVED AFTER SUDIM1S)
!              --------------------------

CALL SULUN1S

!*       7.    INITIALIZE PHYSICS.
!              -------------------

CALL SUPHEC(NULOUT)


!     -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU0YOM1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SU0YOM1S
