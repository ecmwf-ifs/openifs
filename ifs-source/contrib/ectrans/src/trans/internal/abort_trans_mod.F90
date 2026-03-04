! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE ABORT_TRANS_MOD
CONTAINS
SUBROUTINE ABORT_TRANS(CDTEXT)

USE TPM_GEN  , ONLY : NOUT,NERR
USE TPM_DISTR, ONLY : NPROC,MYPROC
USE MPL_MODULE, ONLY : MPL_ABORT
USE SDL_MOD, ONLY : SDL_TRACEBACK, SDL_SRLABORT

IMPLICIT NONE


CHARACTER(LEN=*),INTENT(IN) :: CDTEXT

WRITE(NOUT,'(1X,A)') 'ABORT_TRANS CALLED'

WRITE(NOUT,'(1X,A)') CDTEXT
WRITE(NERR,'(1X,A,1X,I3,1X,A)') 'ABORT! ',MYPROC,CDTEXT
CLOSE(NOUT)
IF (NPROC > 1) THEN
  CALL MPL_ABORT(CDTEXT)
ELSE
  CALL SDL_TRACEBACK
  CALL FLUSH(0)
  CALL SDL_SRLABORT
ENDIF

END SUBROUTINE ABORT_TRANS
END MODULE ABORT_TRANS_MOD
