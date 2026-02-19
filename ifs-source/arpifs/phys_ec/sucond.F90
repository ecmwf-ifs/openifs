! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE SUCOND(YDECND,KULOUT)

!**** *SUCOND* - INITIALIZE COMMON YOECND CONTROLLING LARGE-SCALE COND.

!     PURPOSE.
!     --------
!           INITIALIZE YOECND, THE COMMON THAT CONTROLS THE
!           LARGE-SCALE CONDENSATION ROUTINE OF THE MODEL.

!**   INTERFACE.
!     ----------
!        CALL *SUCOND* FROM *SUPHEC*

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECND

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "IN CORE MODEL"

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-12-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Bechtold    11-Dec-2005 Change value of REPQMI
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOECND   , ONLY : TECND

IMPLICIT NONE

TYPE(TECND)       ,INTENT(INOUT) :: YDECND
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
!      ----------------------------------------------------------------

!*       0.    ARRAY DIMENSIONS
!              ----------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

IF (LHOOK) CALL DR_HOOK('SUCOND',0,ZHOOK_HANDLE)
ASSOCIATE(REPFLM=>YDECND%REPFLM, REPQMI=>YDECND%REPQMI)
REPFLM=1.E-24_JPRB

!REPQMI=1.E-12_JPRB
REPQMI=1.E-8_JPRB

!     -----------------------------------------------------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOECND '')')

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCOND',1,ZHOOK_HANDLE)
END SUBROUTINE SUCOND
