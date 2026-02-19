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

SUBROUTINE SUSCM(KULOUT)

!**** *SUSCM*   - Initialize MODULE YOMSCM controlling extraction of Single-Column profiles.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUSCM(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        MODULE YOMSCM

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        J.M. Piriou.

!     Modifications.
!     --------------
!        Original : 2002-03-17
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMSCM  , ONLY : &
 & LGSCM, NFRSCM, NSCMTS, NSCM_SPACE_S &
 & , GSCM_LON1, GSCM_LON2, GSCM_LAT1, GSCM_LAT2 &
 & , GSCM_RADIUS, NSCM_ADD_SAMPL  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM) ::J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

#include "namscm.nam.h"

!-------------------------------------------------
! 1. Default values.
!-------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSCM',0,ZHOOK_HANDLE)
LGSCM=.FALSE.
NFRSCM=1
NSCMTS=0
NSCM_SPACE_S=0
NSCM_ADD_SAMPL=1
GSCM_LON1=0.0_JPRB
GSCM_LON2=0.0_JPRB
GSCM_LAT1=0.0_JPRB
GSCM_LAT2=0.0_JPRB
GSCM_RADIUS=2.E-04_JPRB
!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

CALL POSNAM(NULNAM,'NAMSCM')
READ(NULNAM,NAMSCM)
!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='(A)') ' MODULE YOMSCM'
WRITE(UNIT=KULOUT,FMT='(A,L5)') '   LGSCM = ',LGSCM
IF(LGSCM) THEN
  WRITE(UNIT=KULOUT,FMT='(3(A,I5),5(A,G12.6))') &
   & '   NFRSCM = ',NFRSCM &
   & ,'NSCM_SPACE_S = ',NSCM_SPACE_S &
   & ,'NSCM_ADD_SAMPL = ',NSCM_ADD_SAMPL &
   & ,'GSCM_LON1 = ',GSCM_LON1 &
   & ,'GSCM_LON2 = ',GSCM_LON2 &
   & ,'GSCM_LAT1 = ',GSCM_LAT1 &
   & ,'GSCM_LAT2 = ',GSCM_LAT2 &
   & ,'GSCM_RADIUS = ',GSCM_RADIUS  
  WRITE(KULOUT,*) '   NSCMTS =  ',NSCMTS(0),(NSCMTS(J),J=1,ABS(NSCMTS(0)))
ENDIF
IF (LHOOK) CALL DR_HOOK('SUSCM',1,ZHOOK_HANDLE)
END SUBROUTINE SUSCM
