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

SUBROUTINE SIPC_READ_MODEL(KINDEX,KDIM,CDJOBNAM_R,KINFOS,PINP_FIELD)  
!****
!**** *SIPC_Read_Model*  

!     Purpose:
!     -------
!     READ model input field and infos from field-specific pool

!**   Interface:
!     ---------
!       *CALL*  *SIPC_Read_Model(...)

!     Input:
!     -----
!          kindex : index of field to be read in total number of input fields
!          kdim   : dimension of field to be read 

!     Output:
!     ------
!          cdjobnam_r   :  Experiment name (character)
!          kinfos(1)    :  Initial date (integer)
!          kinfos(2)    :  Iteration number (integer)
!          kinfos(3)    :  Time since start (integer)
!          pinp_field   :  Field read

!     Workspace:
!     ---------
!     None

!     External
!     ---------
!     SVIPC_read

!     Author:
!     -------
!      S. Valcke      97/08/25

!     Modifications:
!     --------------
! --------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOM_OAS   , ONLY : CPMODINF
USE PAR_SIPC  , ONLY : JPBYTEINT, JPBYTEREA, JPBYTECHA
USE YOM_INPC  , ONLY : MPOOLREAD
USE YOMLUN    , ONLY : NULOUT

! --------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)    :: KDIM 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KINDEX 
CHARACTER (LEN = 3),INTENT(OUT)   :: CDJOBNAM_R
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KINFOS(3)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PINP_FIELD(KDIM)

! --------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISIZEINP, IMRC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!* ---------------------------- Poema verses ----------------------------

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*    1. Initializations
!        ---------------

IF (LHOOK) CALL DR_HOOK('SIPC_READ_MODEL',0,ZHOOK_HANDLE)
WRITE (NULOUT,*) ' '
WRITE (NULOUT,*) &
 & '           ROUTINE SIPC_Read_Model '  
WRITE (NULOUT,*) &
 & '           ***********************'  
WRITE (NULOUT,*) ' '

!*    2. Read encapsulated infos in field-specific shared memory pool
!        (experiment name, initial date, time since start, iteration number)
!        -------------------------------------------------------------------

IF(CPMODINF == 'YES') THEN
  IMRC = 0
  ISIZEINP = 3*JPBYTECHA
  CALL SVIPC_READ(MPOOLREAD(KINDEX), CDJOBNAM_R, ISIZEINP, IMRC)
  ISIZEINP = 3*JPBYTEINT
  CALL SVIPC_READ(MPOOLREAD(KINDEX), KINFOS, ISIZEINP, IMRC)

!*        Find error if any

  IF (IMRC < 0) THEN
    WRITE(NULOUT,*)&
     & 'Problem in reading encapsulated infos for field',&
     & kindex   
    CALL ABOR1('STOP in SIPC_Read_Model')
  ENDIF
ENDIF

!*    3. Read field in field-specific shared memory pool
!        -----------------------------------------------
ISIZEINP = KDIM*JPBYTEREA
CALL SVIPC_READ(MPOOLREAD(KINDEX), PINP_FIELD, ISIZEINP, IMRC)

!*    Find error and STOP if any

IF (IMRC < 0) THEN
  WRITE(NULOUT,*) &
   & 'Problem in reading in the SHM pool for field',&
   & kindex   
  CALL ABOR1('STOP in SIPC_Read_Model')
ENDIF

WRITE(NULOUT,*) 'Read the SHM pool for field:', kindex, cdjobnam_r
WRITE(NULOUT,*)  PINP_FIELD(1), PINP_FIELD(KDIM)

!*    4. End of routine
!        --------------

WRITE (NULOUT,*) ' '
WRITE (NULOUT,*) ' --- End of routine SIPC_Read_model ---'
CALL FLUSH(NULOUT)
IF (LHOOK) CALL DR_HOOK('SIPC_READ_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE SIPC_READ_MODEL
