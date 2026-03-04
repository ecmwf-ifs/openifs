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

SUBROUTINE SIPC_WRITE_MODEL(KINDEX,KDIM,KINFOS,POUT_FIELD)
!****
!**** *SIPC_Write_Model*  

!     Purpose:
!     -------
!     Write model output field and infos to field-specific pool

!**   Interface:
!     ---------
!       *CALL*  *SIPC_Write_Model(...)

!     Input:
!     -----
!        kindex : index of field to be written in total number of output fields
!        kdim   : dimension of field to be written 

!     Output:
!     ------
!        kinfos(1)    :  Initial date (integer)
!        kinfos(2)    :  Iteration number (integer)
!        kinfos(3)    :  Time since start (integer)
!        pout_field   :  Field to be written

!     Workspace:
!     ---------
!     None

!     Externals:
!     ---------
!     SVIPC_write

!     Author:
!     -------
!      S. Valcke      97/08/25

!     Modifications:
!     --------------
! --------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOM_OAS   , ONLY : CPMODINF, CPJOBNAM
USE PAR_SIPC  , ONLY : JPBYTEINT, JPBYTEREA, JPBYTECHA
USE YOM_INPC  , ONLY : MPOOLWRIT
USE YOMLUN    , ONLY : NULOUT

! --------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDEX 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINFOS(3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT_FIELD(KDIM)

! --------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISIZEOUT, IMRC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!* ---------------------------- Poema verses ----------------------------

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*    1. Initializations
!        ---------------

IF (LHOOK) CALL DR_HOOK('SIPC_WRITE_MODEL',0,ZHOOK_HANDLE)
WRITE (NULOUT,*) ' '
WRITE (NULOUT,*) &
 & '           ROUTINE SIPC_Write_Model '  
WRITE (NULOUT,*) &
 & '           *******************     *******'  
WRITE (NULOUT,*) ' '

!*    2. Write encapsulated infos in field-specific shared memory pool
!        (experiment name, initial date, time since start, iteration number)
!        -------------------------------------------------------------------

IF(CPMODINF == 'YES') THEN
  IMRC = 0
  ISIZEOUT = 3*JPBYTECHA
  CALL SVIPC_WRITE(MPOOLWRIT(KINDEX), CPJOBNAM, ISIZEOUT, IMRC)
  ISIZEOUT = 3*JPBYTEINT
  CALL SVIPC_WRITE(MPOOLWRIT(KINDEX), KINFOS, ISIZEOUT, IMRC)

!*        Find error if any

  IF (IMRC < 0) THEN
    WRITE(NULOUT,*)&
     & 'Problem in writing encapsulated infos for field',&
     & kindex   
    CALL ABOR1('STOP in SIPC_Write_Model')
  ENDIF

  WRITE(NULOUT,*) 'Write encapsulated infos for field:',kindex

ENDIF

!*    3. Write field in field-specific shared memory pool
!        -----------------------------------------------
ISIZEOUT = KDIM*JPBYTEREA
CALL SVIPC_WRITE(MPOOLWRIT(KINDEX), POUT_FIELD, ISIZEOUT, IMRC)

!*    Find error and STOP if any

IF (IMRC < 0) THEN
  WRITE(NULOUT,*) &
   & 'Problem in writing in the SHM pool for field',&
   & kindex   
  CALL ABOR1('STOP in SIPC_Write_Model')
ENDIF

WRITE(NULOUT,*) 'Write in the SHM pool for field :',kindex
WRITE(NULOUT,*)  POUT_FIELD(1), POUT_FIELD(KDIM)

!*    4. End of routine
!        --------------

WRITE (NULOUT,*) ' '
WRITE (NULOUT,*)' --- End of routine SIPC_Write_model ---'
CALL FLUSH(NULOUT)
IF (LHOOK) CALL DR_HOOK('SIPC_WRITE_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE SIPC_WRITE_MODEL
