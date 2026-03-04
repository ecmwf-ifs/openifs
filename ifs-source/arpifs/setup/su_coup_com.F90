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

SUBROUTINE SU_COUP_COM(YDCOM,KULOUT)

!**** *SU_COUP_COM* * - ROUTINE TO INITIALIZE VARIABLES FOR THE COM COUPLING

!     PURPOSE.
!     --------
!        SET DEFAULT VALUES, THEN READS NAMELIST NAMCOM

!**   INTERFACE.
!     ----------
!        *CALL* *SU_COUP_COM(KULOUT)*

!     EXPLICIT ARGUMENTS :  KULOUT
!     --------------------

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMCOM

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        Jean-Philippe Piedelievre CNRM/GMGEC/EAC

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 99-02-24
!        Modified : 02-01-11 P. Marquet - Modifications for Climat-V4
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Jan 2010): remove useless variables.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCOM, ONLY : TCOM
USE YOMLUN, ONLY : NULNAM

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCOM),TARGET ,INTENT(INOUT) :: YDCOM
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM), POINTER :: NVCOM
LOGICAL, POINTER :: LOMLDTH

#include "namcom.nam.h"
#include "posnam.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU_COUP_COM',0,ZHOOK_HANDLE)

! Associate pointers for variables in namelist
NVCOM   => YDCOM%NVCOM
LOMLDTH => YDCOM%LOMLDTH

!-----------------------------------------------------------------------

!*     1.   SET DEFAULT VALUES.
!           ------------------

LOMLDTH=.FALSE.
NVCOM=2

!*     2.   MODIFIES DEFAULT VALUES.
!           ------------------------

CALL POSNAM(NULNAM,'NAMCOM')
READ       (NULNAM, NAMCOM)

!*     3.   PRINTS FINAL VALUES.
!           ------------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMCOM '')')
WRITE(UNIT=KULOUT,FMT='('' NVCOM  = '',I5,5X)') NVCOM

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_COUP_COM',1,ZHOOK_HANDLE)
END SUBROUTINE SU_COUP_COM
