! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_NPROC_MOD
!**** MPL_NPROC - return Number of processes 

!        Input optional arguments :
!        -------------------------
!           KCOMM    -  Communicator number if different from MPI_COMM_WORLD 

!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE
PRIVATE
PUBLIC MPL_NPROC

CONTAINS 
FUNCTION MPL_NPROC(KCOMM)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_COMM_SIZE => MPI_COMM_SIZE8
#endif

INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KCOMM
INTEGER(KIND=JPIM) :: MPL_NPROC

INTEGER(KIND=JPIM) :: IERROR,IPROC
LOGICAL   :: LLABORT=.TRUE.

IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE( &
  & CDMESSAGE='MPL_MYRANK: MPL NOT INITIALISED ',LDABORT=LLABORT) 
IF(PRESENT(KCOMM)) THEN
  CALL MPI_COMM_SIZE(KCOMM,IPROC,IERROR)
  MPL_NPROC = IPROC
ELSE
  MPL_NPROC = MPL_NUMPROC 
ENDIF

  
END FUNCTION MPL_NPROC
END MODULE MPL_NPROC_MOD

